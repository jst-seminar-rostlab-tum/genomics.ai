import React from 'react';
import {
  Chip, Stack, TextField, MenuItem,
} from '@mui/material';

function TeamHeaderOptions({
  team, isAdmin, institution, availableInstitutions, setInstitution,
}) {
  const handleInstitutionChange = (event) => {
    setInstitution(event.target.value);
  };

  return (
    <Stack
      direction="row"
      alignItems="flex-end"
      spacing={2}
      sx={{ paddingLeft: '20px' }}
    >
      {(!isAdmin && institution) && <h4>{institution.name}</h4>}
      {!isAdmin && <Chip label={team.visibility} color="primary" />}

      {isAdmin
        && (
          <TextField
            id="select-institution"
            select
            label="Institution"
            value={institution}
            onChange={handleInstitutionChange}
            variant="standard"
          >
            {availableInstitutions.map((institutionOption) => (
              <MenuItem
                key={institutionOption.id}
                value={institutionOption}
              >
                {institutionOption.name}
              </MenuItem>
            ))}
          </TextField>
        )}

    </Stack>
  );
}

export default TeamHeaderOptions;
