import React from 'react';
import { Modal } from 'components/Modal';
import { LearnMoreAtlasComponent } from 'views/Explore/LearnMoreAtlas';
import { Container } from '@mui/material';

/**
 * Modal displaying info about the atlas with the given atlas id
 * @param id Atlas id
 * @param open True if the info modal should be shown
 * @param open Function accepting true or false as a parameter
 */
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
