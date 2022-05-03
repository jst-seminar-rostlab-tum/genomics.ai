import React from 'react';
import { Chip, Stack } from '@mui/material';

function TeamHeaderOptions({
  team, isAdmin, institution,
}) {
  return (
    <Stack
      direction="row"
      alignItems="flex-end"
      spacing={2}
      sx={{ paddingLeft: '20px', paddingBottom: '5px', alignSelf: 'flex-end' }}
    >
      {institution && <h4>{institution.name}</h4>}
      {!isAdmin && <Chip label={team.visibility} color="primary" />}

    </Stack>
  );
}

export default TeamHeaderOptions;
