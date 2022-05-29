/* eslint-disable no-return-assign */
import React, { useState } from 'react';
import Button from 'components/CustomButton';
import TextField from '@mui/material/TextField';
import { Modal, ModalTitle } from 'components/Modal';
import DialogActions from '@mui/material/DialogActions';
import DialogContent from '@mui/material/DialogContent';
import Alert from '@mui/material/Alert';
import InstitutionService from 'shared/services/Institution.service';
import { CircularProgress } from '@mui/material';

export default function InstitutionCreationDialog({ open, handleClose, onCreated }) {
  const [name, setName] = useState('');
  const [country, setCountry] = useState('');
  const [loading, setLoading] = useState(false);
  const [nameError, setNameError] = useState('');
  const [countryError, setCountryError] = useState('');
  const [creationError, setCreationError] = useState('');

  async function create() {
    if (!name) {
      setNameError('Please enter a name.');
    }
    if (!country) {
      setCountryError('Please enter a description.');
    }
    if (!name || !country) return;
    setLoading(true);
    try {
      const newInstitution = await InstitutionService.createInstitution(name, country);
      onCreated(newInstitution);
      handleClose();
      setName('');
      setCountry('');
    } catch (e) {
      setCreationError(e.response.data);
    } finally {
      setLoading(false);
    }
  }

  return (
    <Modal isOpen={open} setOpen={(o) => !o && handleClose()}>
      <ModalTitle>Create an Institution</ModalTitle>
      <DialogContent>
        <TextField
          autoFocus
          label="Name"
          type="text"
          variant="standard"
          required
          error={!!nameError}
          helperText={nameError}
          onChange={(evt) => { setName(evt.target.value); setNameError(''); }}
        />
        <TextField
          label="Country"
          type="text"
          fullWidth
          multiline
          variant="standard"
          required
          error={!!countryError}
          helperText={countryError}
          onChange={(evt) => { setCountry(evt.target.value); setCountryError(''); }}
        />
        {creationError && (
          <Alert severity="error" sx={{ marginTop: '12px' }}>
            {creationError}
          </Alert>
        )}
      </DialogContent>
      <DialogActions>
        <Button type="tertiary" onClick={handleClose}>Cancel</Button>
        {
          loading ? (
            <CircularProgress />
          ) : (
            <Button type="primary" disabled={!name || !country} onClick={() => create()}>
              Create
            </Button>
          )
        }
      </DialogActions>
    </Modal>
  );
}
