import React from 'react';
import { Modal } from 'components/Modal';
import { LearnMoreAtlasComponent } from 'views/Explore/LearnMoreAtlas';
import { Container } from '@mui/material';

function AtlasInfo({ id, open, setOpen }) {
  return (
    <Modal
      isOpen={open}
      setOpen={setOpen}
    >
      <Container>
        <LearnMoreAtlasComponent id={id} />
      </Container>
    </Modal>
  );
}

export default AtlasInfo;
