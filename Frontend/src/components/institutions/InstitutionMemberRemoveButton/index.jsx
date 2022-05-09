import React, { useState } from 'react';
import Button from 'components/CustomButton';
import Dialog from '@mui/material/Dialog';
import DialogActions from '@mui/material/DialogActions';
import DialogContent from '@mui/material/DialogContent';
import DialogContentText from '@mui/material/DialogContentText';
import DialogTitle from '@mui/material/DialogTitle';

import InstitutionService from 'shared/services/Institution.service';

function InstitutionMemberRemoveButton({ institution, member, onRemoved }) {
  const [dialogOpen, setDialogOpen] = useState(false);

  const handleOpenDialog = () => setDialogOpen(true);
  const handleCloseDialog = () => setDialogOpen(false);

  async function remove() {
    await InstitutionService.removeMemberFromInstitution(institution.id, member.id);
    handleCloseDialog();
    onRemoved(institution, member);
  }

  return (
    <>
      <Button type="critical" onClick={handleOpenDialog}>
        Remove
      </Button>
      <Dialog
        open={dialogOpen}
        onClose={handleCloseDialog}
        aria-labelledby="alert-dialog-title"
        aria-describedby="alert-dialog-description"
      >
        <DialogTitle id="alert-dialog-title">
          Remove Member
        </DialogTitle>
        <DialogContent>
          <DialogContentText id="alert-dialog-description">
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
      </Dialog>
    </>
  );
}

export default InstitutionMemberRemoveButton;
