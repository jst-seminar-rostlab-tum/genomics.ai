import React, { useState } from 'react';
import {
  TextField, Dialog, DialogActions, DialogContent,
  DialogContentText, DialogTitle, Button, Fab,
} from '@mui/material';
import PersonAddOutlinedIcon from '@mui/icons-material/PersonAddOutlined';

function TeamInviteButton({ team }) {
  const [dialogOpen, setDialogOpen] = useState(false);
  console.log(dialogOpen);

  const handleCloseDialog = () => setDialogOpen(false);

  const handleTeamInvite = () => {
    handleCloseDialog();
  };

  return (

    <div>
      <Fab
        onClick={() => {
          setDialogOpen(true);
        }}
        color="primary"
        aria-label="add"
        sx={{
          position: 'fixed',
          bottom: '3%',
          right: '2%',
        }}
      >
        <PersonAddOutlinedIcon />
      </Fab>
      <Dialog
        open={dialogOpen}
        onClose={handleCloseDialog}
        aria-labelledby="alert-dialog-title"
        aria-describedby="alert-dialog-description"
        fullWidth
      >
        <DialogTitle id="alert-dialog-title">
          Invite User
        </DialogTitle>
        <DialogContent>
          <DialogContentText id="alert-dialog-description">
            {`Invite a user to your team "${team.name}"`}
          </DialogContentText>
          <TextField
            autoFocus
            margin="dense"
            id="name"
            label="Email Address"
            type="email"
            fullWidth
            variant="standard"
          />
        </DialogContent>
        <DialogActions>
          <Button onClick={handleCloseDialog}>Cancel</Button>
          <Button onClick={handleTeamInvite} autoFocus>
            Confirm
          </Button>
        </DialogActions>
      </Dialog>
    </div>

  );
}

export default TeamInviteButton;
