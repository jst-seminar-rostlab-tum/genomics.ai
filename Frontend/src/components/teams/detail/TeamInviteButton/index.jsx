import React, { useState } from 'react';
import {
  TextField, DialogActions, DialogContent,
  DialogContentText, Fab, Snackbar, Alert,
} from '@mui/material';
import { Modal, ModalTitle } from 'components/Modal';
import Button from 'components/CustomButton';
import PersonAddOutlinedIcon from '@mui/icons-material/PersonAddOutlined';

import TeamService from 'shared/services/Team.service';

function TeamInviteButton({ team }) {
  const [dialogOpen, setDialogOpen] = useState(false);
  const [open, setOpen] = useState(false);
  const [invitedMailAdress, setInvitedMailAdress] = useState('');

  const [mailError, setMailError] = useState('');

  const handleCloseDialog = () => setDialogOpen(false);

  function validateEmail(email) {
    let re = /^[a-zA-Z0-9.!#$%&'*+/=?^_`{|}~-]+@[a-zA-Z0-9-]+(?:\.[a-zA-Z0-9-]+)*$/;
    return re.test(email);
  }

  async function handleTeamInvite() {
    if (!invitedMailAdress || !validateEmail(invitedMailAdress)) {
      setMailError('Please enter a valid e-mail address.');
      return;
    }
    setMailError('');
    try {
      await TeamService.inviteMemberByEmail(team.id, invitedMailAdress);
      handleCloseDialog();
      setOpen(true);
      handleCloseDialog();
    } catch (e) {
      setMailError(e.message);
    }
  }

  const handleClose = (event, reason) => {
    if (reason === 'clickaway') {
      return;
    }
    setOpen(false);
  };

  return (

    <div>
      <Fab
        onClick={() => {
          setDialogOpen(true);
          setInvitedMailAdress('');
          setMailError('');
        }}
        aria-label="add"
        sx={{
          position: 'fixed',
          bottom: '3%',
          right: '2%',
        }}
        style={{
          backgroundColor: '#5676E4'
        }}
      >
        <PersonAddOutlinedIcon style={{ color: 'white' }} />
      </Fab>
      <Modal
        isOpen={dialogOpen}
        setOpen={(o) => !o && handleCloseDialog()}
        fullWidth
      >
        <ModalTitle>
          Invite User
        </ModalTitle>
        <DialogContent>
          <DialogContentText>
            {`Invite a user to your team "${team.name}"`}
          </DialogContentText>
          <TextField
            autoFocus
            required
            margin="dense"
            id="name"
            label="Email Address"
            type="email"
            fullWidth
            variant="standard"
            error={!!mailError}
            helperText={mailError}
            onChange={(evt) => { setInvitedMailAdress(evt.target.value); setMailError(''); }}
          />
        </DialogContent>
        <DialogActions>
          <Button type="tertiary" onClick={handleCloseDialog}>Cancel</Button>
          <Button onClick={handleTeamInvite} autoFocus>
            Confirm
          </Button>
        </DialogActions>
      </Modal>
      <Snackbar
        open={open}
        autoHideDuration={3000}
        onClose={handleClose}
        anchorOrigin={{
          vertical: 'bottom',
          horizontal: 'center',
        }}
      >
        <Alert onClose={handleClose} severity="success" sx={{ width: '100%' }}>
          {'Invite to '}
          {invitedMailAdress}
          {' sent successfully!'}
        </Alert>
      </Snackbar>
    </div>

  );
}

export default TeamInviteButton;
