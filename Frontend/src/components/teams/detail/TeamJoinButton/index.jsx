import React, { useState } from 'react';
import { Modal, ModalTitle } from 'components/Modal';
import DialogActions from '@mui/material/DialogActions';
import DialogContent from '@mui/material/DialogContent';
import DialogContentText from '@mui/material/DialogContentText';
import Button from 'components/CustomButton';
import TeamService from 'shared/services/Team.service';

function TeamJoinButton({ team, onJoin, isDisabled }) {
  const [dialogOpen, setDialogOpen] = useState(false);
  const [errorMessage, setErrorMessage] = useState('');

  const handleOpenDialog = () => setDialogOpen(true);
  const handleCloseDialog = () => setDialogOpen(false);

  async function join() {
    setErrorMessage('');
    try {
      await TeamService.joinTeam(team.id);
      handleCloseDialog();
      onJoin(team);
    } catch (e) {
      setErrorMessage(e.message);
    }
  }

  return (
    <>
      <Button type="primary" disabled={isDisabled} onClick={handleOpenDialog}>Join</Button>
      <Modal
        isOpen={dialogOpen}
        setOpen={(o) => !o && handleCloseDialog()}
      >
        <ModalTitle>
          Join
        </ModalTitle>
        <DialogContent>
          <DialogContentText>
            Do you really want to join the team &quot;
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
          <Button onClick={handleCloseDialog} type="critical">Cancel</Button>
          <Button onClick={() => join()} type="primary" autoFocus>
            Join
          </Button>
        </DialogActions>
      </Modal>
    </>
  );
}

export default TeamJoinButton;
