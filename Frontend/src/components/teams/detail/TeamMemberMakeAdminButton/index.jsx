import React, { useState } from 'react';
import Button from '@mui/material/Button';
import Dialog from '@mui/material/Dialog';
import DialogActions from '@mui/material/DialogActions';
import DialogContent from '@mui/material/DialogContent';
import DialogContentText from '@mui/material/DialogContentText';
import DialogTitle from '@mui/material/DialogTitle';

function TeamMemberMakeAdminButton({ team, member, onMakeAdmin, onRemoveAdmin }) {
  const [dialogOpen, setDialogOpen] = useState(false);

  const handleOpenDialog = () => setDialogOpen(true);
  const handleCloseDialog = () => setDialogOpen(false);

  async function addAdmin() {
    onMakeAdmin(team, member);
    handleCloseDialog();
  }

  async function removeAdmin() {
    onRemoveAdmin(team, member);
    handleCloseDialog();
  }

  const isAlreadyAdmin = team.adminIds ? team.adminIds.includes(member.id) : false;

  return (
    <>
      <Button variant="outlined" onClick={handleOpenDialog} sx={{ marginRight: '20px', width: '145px' }}>
        {isAlreadyAdmin ? 'Remove Admin' : 'Make Admin'}
      </Button>
      <Dialog
        open={dialogOpen}
        onClose={handleCloseDialog}
        aria-labelledby="alert-dialog-title"
        aria-describedby="alert-dialog-description"
      >
        <DialogTitle id="alert-dialog-title">
          {isAlreadyAdmin ? 'Remove Admin' : 'Make Admin'}
        </DialogTitle>
        <DialogContent>
          <DialogContentText id="alert-dialog-description">
            {`Do you really want to ${isAlreadyAdmin ? 'remove' : 'make'} ${member.firstName} 
            ${member.lastName} ${isAlreadyAdmin ? 'as' : ''} an admin of ${team.name}?`}
          </DialogContentText>
        </DialogContent>
        <DialogActions>
          <Button onClick={handleCloseDialog}>Cancel</Button>
          <Button onClick={isAlreadyAdmin ? () => removeAdmin() : () => addAdmin()} color="error" autoFocus>
            {isAlreadyAdmin ? 'Remove' : 'Make'}
          </Button>
        </DialogActions>
      </Dialog>
    </>
  );
}

export default TeamMemberMakeAdminButton;
