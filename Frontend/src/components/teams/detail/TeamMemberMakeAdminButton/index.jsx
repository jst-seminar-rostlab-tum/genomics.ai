import React, { useState } from 'react';
import Button from 'components/CustomButton';
import { Modal, ModalTitle } from 'components/Modal';
import DialogActions from '@mui/material/DialogActions';
import DialogContent from '@mui/material/DialogContent';
import DialogContentText from '@mui/material/DialogContentText';

function TeamMemberMakeAdminButton({
  team, member, onMakeAdmin, onRemoveAdmin,
}) {
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
      <Button type="secondary" onClick={handleOpenDialog}>
        {isAlreadyAdmin ? 'Remove Admin' : 'Make Admin'}
      </Button>
      <Modal
        isOpen={dialogOpen}
        setOpen={(o) => !o && handleCloseDialog()}
        aria-labelledby="alert-dialog-title"
        aria-describedby="alert-dialog-description"
      >
        <ModalTitle id="alert-dialog-title">
          {isAlreadyAdmin ? 'Remove Admin' : 'Make Admin'}
        </ModalTitle>
        <DialogContent>
          <DialogContentText id="alert-dialog-description">
            {`Do you really want to ${isAlreadyAdmin ? 'remove' : 'make'} ${member.firstName} 
            ${member.lastName} ${isAlreadyAdmin ? 'as' : ''} an admin of ${team.name}?`}
          </DialogContentText>
        </DialogContent>
        <DialogActions>
          <Button type="tertiary" onClick={handleCloseDialog}>Cancel</Button>
          <Button type="critical" onClick={isAlreadyAdmin ? () => removeAdmin() : () => addAdmin()} autoFocus>
            {isAlreadyAdmin ? 'Remove' : 'Make'}
          </Button>
        </DialogActions>
      </Modal>
    </>
  );
}

export default TeamMemberMakeAdminButton;
