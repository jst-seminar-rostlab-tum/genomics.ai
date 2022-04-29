import React from "react";
import LabeledSelect from "../LabeledSelect";

// General filter that is needed in all categories
const GeneralFilter = ({  sortBy, onChange }) => {
  const sortItems = [
    { label: "Name", value: "name" },
    { label: "Last updated", value: "lastUpdated" },
  ];

  return (
    <LabeledSelect
      value={sortBy}
      defaultValue={"name"}
      onChange={(event) => onChange("sortBy", event.target.value)}
      items={sortItems}
    />
  );
};

export default GeneralFilter;
