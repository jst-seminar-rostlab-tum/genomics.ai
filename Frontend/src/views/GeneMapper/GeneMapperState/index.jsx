import {
  Box, Container, Step, StepLabel, Stepper,
} from '@mui/material';
import { useState } from 'react';
import AtlasModelChoice from '../AtlasModelChoice/AtlasModelChoice';
import UploadFilePage from '../UploadFilePage';
import { useEffect } from 'react';

function GeneMapperState({ path }) {
  const [selectedAtlas, setSelectedAtlas] = useState('');
  const [selectedModel, setSelectedModel] = useState('');
  const [activeStep, setActiveStep] = useState(0);
  const steps = ['Pick Atlas and Model', 'Choose File and Project details'];

  useEffect(() => {
    setSelectedModel('');
  }, [selectedAtlas]);

  return (
    <Container>
      <Box width="500px" margin="auto" sx={{ marginTop: '1%', marginBottom: '1%' }}>
        <Stepper activeStep={activeStep}>
          {steps.map((labelText, index) => (
            <Step index={index}>
              <StepLabel>
                {labelText}
              </StepLabel>
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
              setSelectedAtlas={setSelectedAtlas}
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
