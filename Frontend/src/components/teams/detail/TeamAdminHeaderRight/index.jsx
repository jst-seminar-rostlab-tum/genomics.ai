import React from 'react';
import { TextField, MenuItem } from '@mui/material';

function TeamAdminHeaderRight({ team, updateVisibility }) {

  return (
    <TextField
      id="select-visibility"
      select
      label="Visibility"
      value={team.visibility}
      onChange={(e) => updateVisibility(e.target.value)}
      variant="standard"
      sx={{ width: '120px' }}
    >
      <MenuItem value="PUBLIC">public</MenuItem>
      <MenuItem value="PRIVATE">private</MenuItem>
      <MenuItem value="BY_INSTITUTION">by institution</MenuItem>
    </TextField>
  );
}

export default TeamAdminHeaderRight;
