import React from 'react';

import {
  Select, MenuItem, Stack, Typography,
} from '@mui/material';

// Probably tem,prorary component to display filter options
function LabeledSelect({
  label, value, onChange, defaultValue, items,
}) {
  return (
    <Stack
      direction="row"
      spacing={2}
      alignItems="center"
      justifyContent="flex-end"
    >
      <Typography sx={{ fontWeight: 'bold' }}>{label}</Typography>
      <Select value={value || defaultValue} onChange={onChange} displayEmpty>
        {items.map((item) => (
          <MenuItem value={item.value}>{item.label}</MenuItem>
        ))}
      </Select>
    </Stack>
  );
}

export default LabeledSelect;
