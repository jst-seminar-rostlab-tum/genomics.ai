import React, { useState } from 'react';
import Button from '@mui/material/Button';
import { Modal, ModalTitle } from 'components/Modal';
import DialogActions from '@mui/material/DialogActions';
import DialogContent from '@mui/material/DialogContent';
import DialogContentText from '@mui/material/DialogContentText';
import CustomButton from 'components/CustomButton';

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
      <CustomButton type="primary" disabled={isDisabled} onClick={handleOpenDialog}>Join</CustomButton>
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
        </DialogContent>
        <DialogActions>
          <Button onClick={handleCloseDialog} color="error">Cancel</Button>
          <Button onClick={() => leave()} autoFocus>
            Join
          </Button>
        </DialogActions>
      </Modal>
    </>
  );
}

export default TeamLeaveButton;
