import React, { useState } from 'react';
import Button from 'components/CustomButton';
import { Modal, ModalTitle } from 'components/Modal';
import DialogActions from '@mui/material/DialogActions';
import DialogContent from '@mui/material/DialogContent';
import DialogContentText from '@mui/material/DialogContentText';

import InstitutionService from 'shared/services/Institution.service';

function InstitutionMemberRemoveButton({ institution, member, updateInstitution }) {
  const [dialogOpen, setDialogOpen] = useState(false);

  const handleOpenDialog = () => setDialogOpen(true);
  const handleCloseDialog = () => setDialogOpen(false);

  async function remove() {
    await InstitutionService.removeMemberFromInstitution(institution.id, member.id);
    updateInstitution();
    handleCloseDialog();
  }

  return (
    <>
      <Button type="critical" onClick={handleOpenDialog}>
        Remove
      </Button>
      <Modal
        isOpen={dialogOpen}
        setOpen={(o) => !o && handleCloseDialog()}
      >
        <ModalTitle>
          Remove Member
        </ModalTitle>
        <DialogContent>
          <DialogContentText>
            Do you really want to remove the member &quot;
            {`${member.firstName} ${member.lastName}`}
            &quot; from the institution &quot;
            {institution.name}
            &quot;?
          </DialogContentText>
        </DialogContent>
        <DialogActions>
          <Button type="tertiary" onClick={handleCloseDialog}>Cancel</Button>
          <Button type="critical" onClick={() => remove()} autoFocus>
            Remove
          </Button>
        </DialogActions>
      </Modal>
    </>
  );
}

export default InstitutionMemberRemoveButton;
