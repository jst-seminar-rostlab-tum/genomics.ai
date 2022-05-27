import React from 'react';
import { TextField, MenuItem } from '@mui/material';

function TeamAdminHeaderRight({ team, setTeam }) {
  const updateVisibility = (event) => {
    setTeam({
      ...team,
      visibility: event.target.value,
    });
  };

  return (
    <TextField
      id="select-visibility"
      select
      label="Visibility"
      value={team.visibility}
      onChange={updateVisibility}
      variant="standard"
      sx={{ width: '120px' }}
    >
      <MenuItem value="public">public</MenuItem>
      <MenuItem value="private">private</MenuItem>
      <MenuItem value="by institution">by institution</MenuItem>
    </TextField>
  );
}

export default TeamAdminHeaderRight;
