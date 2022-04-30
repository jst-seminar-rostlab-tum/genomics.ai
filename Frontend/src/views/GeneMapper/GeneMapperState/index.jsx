import {
  Box, Container, Step, StepLabel, Stepper,
} from '@mui/material';
import { useState } from 'react';
import AtlasModelChoice from '../AtlasModelChoice/AtlasModelChoice';
import UploadFilePage from '../UploadFilePage';

function GeneMapperState({ path }) {
  const [selectedAtlas, setSelectedAtlas] = useState('');
  const [selectedModel, setSelectedModel] = useState('');
  const [activeStep, setActiveStep] = useState(0);
  const steps = ['Pick Atlas and Model', 'Choose File and Project details'];

  return (
    <Container>
      <Box width="500px" margin="auto" sx={{ paddingTop: '4%', marginBottom: '3%' }}>
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
