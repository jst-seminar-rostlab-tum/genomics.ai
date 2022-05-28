import React, { useState, useEffect } from 'react';
import Button from 'components/CustomButton';
import { Modal, ModalTitle } from 'components/Modal';
import DialogActions from '@mui/material/DialogActions';
import DialogContent from '@mui/material/DialogContent';
import DialogContentText from '@mui/material/DialogContentText';
import Tooltip from '@mui/material/Tooltip';

import { useAuth } from 'shared/context/authContext';
import TeamService from 'shared/services/Team.service';

function TeamLeaveButton({ team, onLeft }) {
  const [dialogOpen, setDialogOpen] = useState(false);
  const [errorMessage, setErrorMessage] = useState('');
  const [user] = useAuth();

  const handleOpenDialog = () => setDialogOpen(true);
  const handleCloseDialog = () => setDialogOpen(false);

  const [onlyAdmin, setOnlyAdmin] = useState();
  useEffect(() => {
    setOnlyAdmin(team.adminIds.length === 1 && team.adminIds[0] === user._id);
  }, [team, user]);

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
      <Tooltip title={onlyAdmin ? 'You are the only admin of this team.' : ''}>
        <div>
          <Button
            type="critical"
            onClick={handleOpenDialog}
            disabled={onlyAdmin}
          >
            Leave
          </Button>
        </div>
      </Tooltip>
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
