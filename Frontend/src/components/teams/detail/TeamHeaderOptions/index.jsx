import React, { useState } from 'react';
import {
  Chip, Stack, TextField, MenuItem, Button, Dialog, DialogActions, DialogContent,
  DialogContentText, DialogTitle,
} from '@mui/material';

function TeamHeaderOptions({
  team, isAdmin, institution, availableInstitutions, setInstitution,
}) {
  const [dialogOpen, setDialogOpen] = useState(false);
  const [targetInstitution, setTargetInstitution] = useState({});

  const handleCloseDialog = () => setDialogOpen(false);

  const handleOpenDialog = (event) => {
    setDialogOpen(true);
    setTargetInstitution(event.target.value);
  };

  const handleInstitutionChange = () => {
    setInstitution(targetInstitution);
    handleCloseDialog();
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
            onChange={handleOpenDialog}
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

      <Dialog
        open={dialogOpen}
        onClose={handleCloseDialog}
        aria-labelledby="alert-dialog-title"
        aria-describedby="alert-dialog-description"
      >
        <DialogTitle id="alert-dialog-title">
          Change Institution
        </DialogTitle>
        <DialogContent>
          <DialogContentText id="alert-dialog-description">
            {`Do you really want to move the team ${team.name} to the institution ${targetInstitution.name}?`}
          </DialogContentText>
        </DialogContent>
        <DialogActions>
          <Button onClick={handleCloseDialog}>Cancel</Button>
          <Button onClick={handleInstitutionChange} color="error" autoFocus>
            Confirm
          </Button>
        </DialogActions>
      </Dialog>

    </Stack>
  );
}

export default TeamHeaderOptions;
