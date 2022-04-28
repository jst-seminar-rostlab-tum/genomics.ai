import React from "react";

import { Select, MenuItem } from "@mui/material";

const LabeledSelect = ({ value, onChange, defaultValue, items }) => {
  return <Select value={value || defaultValue} onChange={onChange} displayEmpty>
    {items.map(({ label, value }) => (
      <MenuItem value={value}>{label}</MenuItem>
    ))}
  </Select>;
};

export default LabeledSelect;
