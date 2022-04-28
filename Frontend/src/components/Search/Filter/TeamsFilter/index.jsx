import React, { useEffect } from "react";
import LabeledSelect from "../LabeledSelect";

const TeamsFilter = ({ visibility, onChange }) => {
  const visibilityItems = [
    { label: "All", value: "" },
    { label: "Public", value: "public" },
    { label: "Private", value: "private" },
    { label: "By institution", value: "byInstitution" },
  ];

  useEffect(
    () => () => {
      onChange("visibility", "");
    },
    []
  );

  return (
    <LabeledSelect
      value={visibility}
      defaultValue={""}
      onChange={(event) => onChange("visibility", event.target.value)}
      items={visibilityItems}
    />
  );
};

export default TeamsFilter;
