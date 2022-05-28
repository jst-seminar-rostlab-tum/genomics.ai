import React, { useState, useEffect } from 'react';
import Button from 'components/CustomButton';
import { Modal, ModalTitle } from 'components/Modal';
import DialogActions from '@mui/material/DialogActions';
import DialogContent from '@mui/material/DialogContent';
import DialogContentText from '@mui/material/DialogContentText';
import Tooltip from '@mui/material/Tooltip';

import { useAuth } from 'shared/context/authContext';
import InstitutionService from 'shared/services/Institution.service';

function InstitutionLeaveButton({ institution, onLeft }) {
  const [dialogOpen, setDialogOpen] = useState(false);
  const [errorMessage, setErrorMessage] = useState('');
  const [user] = useAuth();

  const handleOpenDialog = () => setDialogOpen(true);
  const handleCloseDialog = () => setDialogOpen(false);

  const [onlyAdmin, setOnlyAdmin] = useState();
  useEffect(() => {
    setOnlyAdmin(institution.adminIds.length === 1 && institution.adminIds[0] === user._id);
  }, [institution, user]);

  async function leave() {
    setErrorMessage('');
    try {
      await InstitutionService.leaveInstitution(institution.id);
      handleCloseDialog();
      onLeft(institution);
    } catch (e) {
      setErrorMessage(e.message);
    }
  }

  return (
    <>
      <Tooltip title={onlyAdmin ? 'You are the only admin of this institution.' : ''}>
        <div>
          <Button
            type="critical"
            onClick={handleOpenDialog}
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
          Leave institution
        </ModalTitle>
        <DialogContent>
          <DialogContentText>
            Do you really want to leave the institution &quot;
            {institution.name}
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
          <Button onClick={handleCloseDialog}>Cancel</Button>
          <Button onClick={() => leave()} type="critical" autoFocus>
            Leave
          </Button>
        </DialogActions>
      </Modal>
    </>
  );
}

export default InstitutionLeaveButton;
