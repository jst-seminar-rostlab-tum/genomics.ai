import React from 'react';
import { Modal } from 'components/Modal';
import { LearnMoreModelComponent } from 'views/References/LearnMoreModel';
import { Container } from '@mui/material';

/**
 * Modal displaying info about the model with the given model id
 * @param id Model id
 * @param open True if the info modal should be shown
 * @param open Function accepting true or false as a parameter
 */
function ModelInfo({ id, open, setOpen }) {
  return (
    <Modal
      isOpen={open}
      setOpen={setOpen}
    >
      <Container>
        <LearnMoreModelComponent id={id} />
      </Container>
    </Modal>
  );
}

export default ModelInfo;
