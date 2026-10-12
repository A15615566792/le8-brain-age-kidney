from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def dates(values: pd.Series) -> pd.Series:
    return pd.to_datetime(values, errors="coerce").dt.normalize()


def second_persistent_date(group: pd.DataFrame, column: str, threshold: float, lower: bool) -> pd.Timestamp:
    values = pd.to_numeric(group[column], errors="coerce")
    eligible = values.lt(threshold) if lower else values.ge(threshold)
    observed = np.sort(group.loc[eligible & group["measurement_date"].notna(), "measurement_date"].unique())
    for first in observed:
        later = observed[observed >= first + np.timedelta64(90, "D")]
        if later.size:
            return pd.Timestamp(later[0])
    return pd.NaT


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("participants", type=Path)
    parser.add_argument("coded_events", type=Path)
    parser.add_argument("laboratory_measurements", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--administrative-censor-date", default="2021-12-31")
    args = parser.parse_args()

    participants = pd.read_csv(args.participants, low_memory=False)
    try:
        events = pd.read_csv(args.coded_events, low_memory=False)
    except pd.errors.EmptyDataError:
        events = pd.DataFrame(columns=["participant_id", "event_date"])
    laboratory = pd.read_csv(args.laboratory_measurements, low_memory=False)
    required_participants = {"participant_id", "index_date"}
    required_events = {"participant_id", "event_date"}
    required_laboratory = {"participant_id", "measurement_date", "egfr", "uacr"}
    if not required_participants.issubset(participants.columns):
        raise ValueError(f"participants requires {sorted(required_participants)}")
    if not required_events.issubset(events.columns):
        raise ValueError(f"coded_events requires {sorted(required_events)}")
    if not required_laboratory.issubset(laboratory.columns):
        raise ValueError(f"laboratory_measurements requires {sorted(required_laboratory)}")

    participants["participant_id"] = participants["participant_id"].astype("string")
    events["participant_id"] = events["participant_id"].astype("string")
    laboratory["participant_id"] = laboratory["participant_id"].astype("string")
    participants["index_date"] = dates(participants["index_date"])
    for column in ["death_date", "loss_date", "coverage_start", "coverage_end"]:
        if column not in participants:
            participants[column] = pd.NaT
        else:
            participants[column] = dates(participants[column])
    events["event_date"] = dates(events["event_date"])
    if "qualifying_ckd" in events:
        events = events[pd.to_numeric(events["qualifying_ckd"], errors="coerce").eq(1)]
    coded = events.dropna(subset=["event_date"]).groupby("participant_id")["event_date"].min()
    coded = coded.astype("datetime64[ns]")

    laboratory["measurement_date"] = dates(laboratory["measurement_date"])
    laboratory = laboratory.sort_values(["participant_id", "measurement_date"])
    lab_dates = []
    for participant_id, group in laboratory.groupby("participant_id", sort=False):
        egfr_date = second_persistent_date(group, "egfr", 60.0, True)
        uacr_date = second_persistent_date(group, "uacr", 30.0, False)
        eligible = [value for value in [egfr_date, uacr_date] if pd.notna(value)]
        lab_dates.append((participant_id, min(eligible) if eligible else pd.NaT))
    laboratory_event = pd.Series(dict(lab_dates), name="laboratory_event_date", dtype="datetime64[ns]")

    result = participants.copy()
    result["coded_event_date"] = result["participant_id"].map(coded)
    result["laboratory_event_date"] = result["participant_id"].map(laboratory_event)
    result["composite_ckd_date"] = result[["coded_event_date", "laboratory_event_date"]].min(axis=1)
    result["prevalent_ckd"] = result["composite_ckd_date"].le(result["index_date"])
    admin = pd.Timestamp(args.administrative_censor_date)
    censor = pd.DataFrame({
        "administrative": pd.Series(admin, index=result.index),
        "death": result["death_date"],
        "loss": result["loss_date"],
        "coverage": result["coverage_end"],
    }).min(axis=1)
    incident_date = result["composite_ckd_date"].where(result["composite_ckd_date"].gt(result["index_date"]))
    result["kidney_event"] = incident_date.notna() & incident_date.le(censor) & ~result["prevalent_ckd"]
    result["kidney_outcome_date"] = incident_date.where(result["kidney_event"])
    result["end_date"] = pd.concat([result["kidney_outcome_date"], censor], axis=1).min(axis=1)
    result["followup_years"] = (result["end_date"] - result["index_date"]).dt.days / 365.25
    result["death_competing"] = result["death_date"].eq(result["end_date"]) & ~result["kidney_event"]
    result["eligible_followup"] = result["index_date"].notna() & result["end_date"].gt(result["index_date"])
    if result["coverage_start"].notna().any():
        coverage_eligible = result["coverage_start"].isna() | result["coverage_start"].le(result["index_date"])
        result["eligible_followup"] &= coverage_eligible
    columns = [
        "participant_id", "prevalent_ckd", "kidney_event", "kidney_outcome_date",
        "end_date", "followup_years", "death_competing", "eligible_followup",
    ]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    result[columns].to_csv(args.output, index=False)


if __name__ == "__main__":
    main()
