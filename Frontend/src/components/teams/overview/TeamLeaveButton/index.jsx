import React, { useState } from 'react';
import Button from 'components/CustomButton';
import { Modal, ModalTitle } from 'components/Modal';
import DialogActions from '@mui/material/DialogActions';
import DialogContent from '@mui/material/DialogContent';
import DialogContentText from '@mui/material/DialogContentText';

import TeamService from 'shared/services/Team.service';

function TeamLeaveButton({ team, onLeft }) {
  const [dialogOpen, setDialogOpen] = useState(false);
  const [errorMessage, setErrorMessage] = useState('');

  const handleOpenDialog = () => setDialogOpen(true);
  const handleCloseDialog = () => setDialogOpen(false);

  async function leave() {
    setErrorMessage('');
    try {
      await TeamService.leaveTeam(team.id);
      handleCloseDialog();
      onLeft(team);
    } catch (e) {
      setErrorMessage(e.message);
    }
  }

  return (
    <>
      <Button
        type="critical"
        onClick={handleOpenDialog}
      >
        Leave
      </Button>
      <Modal
        isOpen={dialogOpen}
        setOpen={(o) => !o && handleCloseDialog()}
      >
        <ModalTitle>
          Leave
        </ModalTitle>
        <DialogContent>
          <DialogContentText>
            Do you really want to leave the team &quot;
            {team.name}
            &quot;?
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
          <Button onClick={handleCloseDialog} type="tertiary">Cancel</Button>
          <Button onClick={() => leave()} type="critical" autoFocus>
            Leave
          </Button>
        </DialogActions>
      </Modal>
    </>
  );
}

export default TeamLeaveButton;
