/* eslint-disable no-return-assign */
import React, { useState } from 'react';
import Button from 'components/CustomButton';
import TextField from '@mui/material/TextField';
import { Modal, ModalTitle } from 'components/Modal';
import DialogActions from '@mui/material/DialogActions';
import DialogContent from '@mui/material/DialogContent';
import TeamService from 'shared/services/Team.service';
import InstitutionChoice from 'components/institutions/InstitutionChoice';
import { CircularProgress } from '@mui/material';
import { useInstitutions } from 'shared/context/institutionContext';

export default function TeamCreationDialog({ open, handleClose, onCreated }) {
  const [name, setName] = useState('');
  const [description, setDescription] = useState('');
  const [loading, setLoading] = useState(false);
  const [nameError, setNameError] = useState('');
  const [descriptionError, setDescriptionError] = useState('');
  const [institutionId, setInstitutionId] = useState(null);
  const { institutions } = useInstitutions();

  async function create() {
    if (!name) {
      setNameError('Please enter a name.');
    }
    if (!description) {
      setDescriptionError('Please enter a description.');
    }
    if (!name || !description) return;
    setLoading(true);
    const newTeam = await TeamService.createTeam(name, description, institutionId);
    onCreated(newTeam);
    setLoading(false);
    handleClose();
  }

  if (institutions.length === 0) {
    return (
      <>
        <b>You are not an admin of any institution.</b>
        <p>
          To create a team, you need to be the admin of an institution.
          That is because you can only create teams within institutions directly.
        </p>
      </>
    );
  }

  return (
    <Modal
      isOpen={open}
      setOpen={(o) => !o && handleClose()}>
      <ModalTitle>Create a Team</ModalTitle>
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
        <br />
        <br />
        <InstitutionChoice onChoiceChange={setInstitutionId} />
      </DialogContent>
      <DialogActions>
        <Button onClick={handleClose} type="tertiary">
          Cancel
        </Button>
        {
          loading ? (
            <CircularProgress />
          ) : (
            <Button onClick={() => create()} disabled={!name || !institutionId}>
              Create
            </Button>
          )
        }
      </DialogActions>
    </Modal >
  );
}
