import React, { useState } from 'react';
import Button from 'components/CustomButton';
import { Modal, ModalTitle } from 'components/Modal';
import DialogActions from '@mui/material/DialogActions';
import DialogContent from '@mui/material/DialogContent';
import DialogContentText from '@mui/material/DialogContentText';

import InstitutionService from 'shared/services/Institution.service';

function InstitutionMemberMakeAdminButton({
  institution, member, updateInstitution,
}) {
  const [dialogOpen, setDialogOpen] = useState(false);

  const handleOpenDialog = () => setDialogOpen(true);
  const handleCloseDialog = () => setDialogOpen(false);

  async function addAdmin() {
    await InstitutionService.makeInstitutionAdmin(institution.id, member.id);
    updateInstitution();
    handleCloseDialog();
  }

  async function removeAdmin() {
    await InstitutionService.removeInstitutionAdmin(institution.id, member.id);
    updateInstitution();
    handleCloseDialog();
  }

  const isAlreadyAdmin = institution.adminIds ? institution.adminIds.includes(member.id) : false;

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
            ${member.lastName} ${isAlreadyAdmin ? 'as' : ''} an admin of ${institution.name}?`}
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

export default InstitutionMemberMakeAdminButton;
