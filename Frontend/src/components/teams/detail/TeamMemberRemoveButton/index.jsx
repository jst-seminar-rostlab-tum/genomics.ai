import React, { useState } from 'react';
import Button from '@mui/material/Button';
import Dialog from '@mui/material/Dialog';
import DialogActions from '@mui/material/DialogActions';
import DialogContent from '@mui/material/DialogContent';
import DialogContentText from '@mui/material/DialogContentText';
import DialogTitle from '@mui/material/DialogTitle';

import TeamService from 'shared/services/Team.service';

function TeamMemberRemoveButton({ team, member, onRemoved }) {
  const [dialogOpen, setDialogOpen] = useState(false);

  const handleOpenDialog = () => setDialogOpen(true);
  const handleCloseDialog = () => setDialogOpen(false);

  async function remove() {
    await TeamService.removeMemberFromTeam(team.id, member.id);
    handleCloseDialog();
    onRemoved(team, member);
  }

  return (
    <>
      <Button variant="outlined" color="error" onClick={handleOpenDialog}>
        Remove
      </Button>
      <Dialog
        open={dialogOpen}
        onClose={handleCloseDialog}
        aria-labelledby="alert-dialog-title"
        aria-describedby="alert-dialog-description"
      >
        <DialogTitle id="alert-dialog-title">
          Remove Member
        </DialogTitle>
        <DialogContent>
          <DialogContentText id="alert-dialog-description">
            Do you really want to remove the member &quot;
            {`${member.firstName} ${member.lastName}`}
            &quot; from the team &quot;
            {team.name}
            &quot;?
          </DialogContentText>
        </DialogContent>
        <DialogActions>
          <Button onClick={handleCloseDialog}>Cancel</Button>
          <Button onClick={() => remove()} color="error" autoFocus>
            Remove
          </Button>
        </DialogActions>
      </Dialog>
    </>
  );
}

export default TeamMemberRemoveButton;
