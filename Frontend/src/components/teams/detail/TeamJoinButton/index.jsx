import React, { useState } from 'react';
import Button from '@mui/material/Button';
import Dialog from '@mui/material/Dialog';
import DialogActions from '@mui/material/DialogActions';
import DialogContent from '@mui/material/DialogContent';
import DialogContentText from '@mui/material/DialogContentText';
import DialogTitle from '@mui/material/DialogTitle';

function TeamLeaveButton({ team, onJoin, isDisabled }) {
  const [dialogOpen, setDialogOpen] = useState(false);

  const handleOpenDialog = () => setDialogOpen(true);
  const handleCloseDialog = () => setDialogOpen(false);

  function leave() {
    handleCloseDialog();
    onJoin(team);
  }

  return (
    <>
      <Button variant="outlined" disabled={isDisabled} onClick={handleOpenDialog}>
        Join
      </Button>
      <Dialog
        open={dialogOpen}
        onClose={handleCloseDialog}
        aria-labelledby="alert-dialog-title"
        aria-describedby="alert-dialog-description"
      >
        <DialogTitle id="alert-dialog-title">
          Join
        </DialogTitle>
        <DialogContent>
          <DialogContentText id="alert-dialog-description">
            Do you really want to join the team &quot;
            {team.name}
            &quot;?
          </DialogContentText>
        </DialogContent>
        <DialogActions>
          <Button onClick={handleCloseDialog} color="critical">Cancel</Button>
          <Button onClick={() => leave()} autoFocus>
            Join
          </Button>
        </DialogActions>
      </Dialog>
    </>
  );
}

export default TeamLeaveButton;
