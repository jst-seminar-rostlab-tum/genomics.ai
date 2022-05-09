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

  const handleOpenDialog = () => setDialogOpen(true);
  const handleCloseDialog = () => setDialogOpen(false);

  async function leave() {
    await InstitutionService.leaveInstitution(institution);
    handleCloseDialog();
    onLeft(institution);
  }

  return (
    <>
      <Button type="critical" onClick={handleOpenDialog}>
        Leave
      </Button>
      <Dialog
        open={dialogOpen}
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
