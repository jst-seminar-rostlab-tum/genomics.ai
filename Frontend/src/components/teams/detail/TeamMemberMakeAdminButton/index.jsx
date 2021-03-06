import React, { useState } from 'react';
import Button from 'components/CustomButton';
import { Modal, ModalTitle } from 'components/Modal';
import DialogActions from '@mui/material/DialogActions';
import DialogContent from '@mui/material/DialogContent';
import DialogContentText from '@mui/material/DialogContentText';

import TeamService from 'shared/services/Team.service';

function TeamMemberMakeAdminButton({
  team, member, updateTeam
}) {
  const [dialogOpen, setDialogOpen] = useState(false);
  const [errorMessage, setErrorMessage] = useState('');

  const handleOpenDialog = () => setDialogOpen(true);
  const handleCloseDialog = () => setDialogOpen(false);

  async function addAdmin() {
    setErrorMessage('');
    try {
      await TeamService.makeTeamAdmin(team.id, member.id);
      handleCloseDialog();
      updateTeam();
    } catch (e) {
      setErrorMessage(e.message);
    }
  }

  async function removeAdmin() {
    setErrorMessage('');
    try {
      await TeamService.removeTeamAdmin(team.id, member.id);
      updateTeam();
      handleCloseDialog();
    } catch (e) {
      setErrorMessage(e.message);
    }
  }

  const isAlreadyAdmin = team.adminIds ? team.adminIds.includes(member.id) : false;

  return (
    <>
      <Button type="primary" onClick={handleOpenDialog}>
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
          {
            errorMessage && (
              <DialogContentText id="alert-dialog-description" color="error">
                {errorMessage}
              </DialogContentText>
            )
          }
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
