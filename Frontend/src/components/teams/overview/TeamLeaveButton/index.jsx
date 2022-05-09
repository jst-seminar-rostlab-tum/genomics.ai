import React, { useState } from 'react';
import Button from 'components/CustomButton';
import Dialog from '@mui/material/Dialog';
import DialogActions from '@mui/material/DialogActions';
import DialogContent from '@mui/material/DialogContent';
import DialogContentText from '@mui/material/DialogContentText';
import DialogTitle from '@mui/material/DialogTitle';

import TeamService from 'shared/services/Team.service';

function TeamLeaveButton({ team, onLeft }) {
  const [dialogOpen, setDialogOpen] = useState(false);

  const handleOpenDialog = () => setDialogOpen(true);
  const handleCloseDialog = () => setDialogOpen(false);

  async function leave() {
    await TeamService.leaveTeam(team.id);
    handleCloseDialog();
    onLeft(team);
  }

  return (
    <>
      <Button variant="outlined" type="critical" onClick={handleOpenDialog}>
        Leave
      </Button>
      <Dialog
        open={dialogOpen}
        onClose={handleCloseDialog}
        aria-labelledby="alert-dialog-title"
        aria-describedby="alert-dialog-description"
      >
        <DialogTitle id="alert-dialog-title">
          Leave
        </DialogTitle>
        <DialogContent>
          <DialogContentText id="alert-dialog-description">
            Do you really want to leave the team &quot;
            {team.name}
            &quot;?
          </DialogContentText>
        </DialogContent>
        <DialogActions>
          <Button onClick={handleCloseDialog}>Cancel</Button>
          <Button onClick={() => leave()} color="error" autoFocus>
            Leave
          </Button>
        </DialogActions>
      </Dialog>
    </>
  );
}

export default TeamLeaveButton;
