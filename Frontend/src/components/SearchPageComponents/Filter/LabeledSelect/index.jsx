import React from "react";

import { Select, MenuItem, Stack, Typography } from "@mui/material";

// Probably tem,prorary component to display filter options
const LabeledSelect = ({ label, value, onChange, defaultValue, items }) => {
  return (
    <Stack
      direction="row"
      spacing={2}
      alignItems="center"
      justifyContent="flex-end"
    >
      <Typography sx={{ fontWeight: "bold" }}>{label}</Typography>
      <Select value={value || defaultValue} onChange={onChange} displayEmpty>
        {items.map(({ label, value }) => (
          <MenuItem value={value}>{label}</MenuItem>
        ))}
      </Select>
    </Stack>
  );
};

export default LabeledSelect;
