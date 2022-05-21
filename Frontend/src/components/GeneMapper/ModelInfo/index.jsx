import React from 'react';
import { Modal } from 'components/Modal';
import { LearnMoreModelComponent } from 'views/Explore/LearnMoreModel';
import { Container } from '@mui/material';

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
