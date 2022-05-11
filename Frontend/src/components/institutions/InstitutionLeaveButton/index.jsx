import React, { useState } from 'react';
import Button from 'components/CustomButton';
import Dialog from '@mui/material/Dialog';
import DialogActions from '@mui/material/DialogActions';
import DialogContent from '@mui/material/DialogContent';
import DialogContentText from '@mui/material/DialogContentText';
import DialogTitle from '@mui/material/DialogTitle';

import InstitutionService from 'shared/services/Institution.service';

function InstitutionLeaveButton({ institution, onLeft }) {
  const [dialogOpen, setDialogOpen] = useState(false);
  const [errorMessage, setErrorMessage] = useState('');

  const handleOpenDialog = () => setDialogOpen(true);
  const handleCloseDialog = () => setDialogOpen(false);

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
      <Button type="critical" onClick={handleOpenDialog}>
        Leave
      </Button>
      <Dialog
        open={dialogOpen}
        onClick={(e) => e.stopPropagation()}
        onClose={handleCloseDialog}
        aria-labelledby="alert-dialog-title"
        aria-describedby="alert-dialog-description"
      >
        <DialogTitle id="alert-dialog-title">
          Leave institution
        </DialogTitle>
        <DialogContent>
          <DialogContentText id="alert-dialog-description">
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
      </Dialog>
    </>
  );
}

export default InstitutionLeaveButton;
