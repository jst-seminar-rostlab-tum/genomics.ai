/* eslint-disable no-return-assign */
import React, { useState } from 'react';
import Button from '@mui/material/Button';
import LoadingButton from '@mui/lab/LoadingButton';
import TextField from '@mui/material/TextField';
import Dialog from '@mui/material/Dialog';
import DialogActions from '@mui/material/DialogActions';
import DialogContent from '@mui/material/DialogContent';
import DialogTitle from '@mui/material/DialogTitle';
import InstitutionService from 'shared/services/Institution.service';

export default function InstitutionCreationDialog({ open, handleClose, onCreated }) {
  const [name, setName] = useState('');
  const [description, setDescription] = useState('');
  const [loading, setLoading] = useState(false);
  const [nameError, setNameError] = useState('');
  const [descriptionError, setDescriptionError] = useState('');

  async function create() {
    if (!name) {
      setNameError('Please enter a name.');
    }
    if (!description) {
      setDescriptionError('Please enter a description.');
    }
    if (!name || !description) return;
    setLoading(true);
    const newInstitution = await InstitutionService.createInstitution(name, description);
    onCreated(newInstitution);
    setLoading(false);
    handleClose();
    setName('');
    setDescription('');
  }

  return (
    <Dialog open={open} onClose={handleClose}>
      <DialogTitle>Create an Institution</DialogTitle>
      <DialogContent>
        <TextField
          autoFocus
          margin="dense"
          id="name"
          label="Name"
          type="text"
          variant="standard"
          required
          error={!!nameError}
          helperText={nameError}
          onChange={(evt) => { setName(evt.target.value); setNameError(''); }}
        />
        <TextField
          margin="dense"
          id="description"
          label="Description"
          type="text"
          fullWidth
          multiline
          maxRows={3}
          variant="standard"
          required
          error={!!descriptionError}
          helperText={descriptionError}
          onChange={(evt) => { setDescription(evt.target.value); setDescriptionError(''); }}
        />
      </DialogContent>
      <DialogActions>
        <Button onClick={handleClose}>Cancel</Button>
        <LoadingButton onClick={() => create()} loading={loading}>Create</LoadingButton>
      </DialogActions>
    </Dialog>
  );
}
