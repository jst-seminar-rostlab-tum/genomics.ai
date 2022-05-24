import {
  Box, Container, Step, StepButton, Stepper,
} from '@mui/material';
import { useState } from 'react';
import AtlasModelChoice from '../AtlasModelChoice/AtlasModelChoice';
import UploadFilePage from '../UploadFilePage';
import { useEffect } from 'react';

function GeneMapperState({ path, defaultSelectedAtlas = '', defaultSelectedModel = '' }) {
  const [selectedAtlas, setSelectedAtlas] = useState(defaultSelectedAtlas);
  const [selectedModel, setSelectedModel] = useState(defaultSelectedModel);
  const [activeStep, setActiveStep] = useState(0);
  const steps = ['Pick Atlas and Model', 'Choose File and Project details'];

  const handleSelectAtlas = (atlas) => {
    setSelectedAtlas(atlas);
    setSelectedModel('');
  }

  const handleStep = (step) => () => {
    setActiveStep(step);
  }

  return (
    <Container>
      <Box width="500px" margin="auto" sx={{ marginTop: '1%', marginBottom: '1%' }}>
        <Stepper activeStep={activeStep}>
          {steps.map((labelText, index) => (
            <Step index={index}>
              <StepButton color="inherit" onClick={handleStep(index)}>
                {labelText}
              </StepButton>
            </Step>
          ))}
        </Stepper>
      </Box>
      {
        activeStep === 0
          ? (
            <AtlasModelChoice
              path={path}
              selectedAtlas={selectedAtlas}
              selectedModel={selectedModel}
              activeStep={activeStep}
              steps={steps}
              setSelectedAtlas={handleSelectAtlas}
              setSelectedModel={setSelectedModel}
              setActiveStep={setActiveStep}
              compatibleModels={selectedAtlas ? selectedAtlas.compatibleModels : []}
            />
          )
          : (
            <UploadFilePage
              path={path}
              selectedAtlas={selectedAtlas}
              selectedModel={selectedModel}
              setActiveStep={setActiveStep}
            />
          )
      }
    </Container>
  );
}

export default GeneMapperState;
