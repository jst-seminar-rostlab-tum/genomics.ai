import React, { useState } from 'react';
import Button from 'components/CustomButton';
import { Modal, ModalTitle } from 'components/Modal';
import DialogActions from '@mui/material/DialogActions';
import DialogContent from '@mui/material/DialogContent';
import DialogContentText from '@mui/material/DialogContentText';

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
