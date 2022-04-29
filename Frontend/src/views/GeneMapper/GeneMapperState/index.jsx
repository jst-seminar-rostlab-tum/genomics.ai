import { Box, Container, Step, StepLabel, Stepper } from '@mui/material';
import { useState } from 'react';
import AtlasModelChoice from '../AtlasModelChoice/AtlasModelChoice';
import UploadFilePage from '../UploadFilePage';

function GeneMapperState(props) {
    const [selectedAtlas, setSelectedAtlas] = useState("");
    const [selectedModel, setSelectedModel] = useState("");
    const [activeStep, setActiveStep] = useState(0);
    const steps = ["Pick Atlas and Model", "Choose File and Project details"];

    const functions = {setSelectedAtlas, setSelectedModel, setActiveStep};

    return (
        <Container>
            <Box width="500px" margin="auto" sx={{ marginTop:'5%', marginBottom:'3%'}}>
                <Stepper activeStep={activeStep}>
                    {steps.map((labelText, index) => {
                        return (
                        <Step index={index}>
                            <StepLabel>
                                {labelText}
                            </StepLabel>
                        </Step>
                        )
                    })}
                </Stepper>
            </Box>
            { 
                activeStep === 0 ? 
                <AtlasModelChoice 
                    selectedAtlas={selectedAtlas}
                    selectedModel={selectedModel}
                    activeStep={activeStep}
                    steps={steps}
                    setSelectedAtlas={setSelectedAtlas}
                    setSelectedModel={setSelectedModel}
                    setActiveStep={setActiveStep}
                /> :
                <UploadFilePage
                    selectedAtlas={selectedAtlas}
                    selectedModel={selectedModel}
                    activeStep={activeStep}
                    steps={steps}
                    setActiveStep={setActiveStep}
                />
            }
        </Container>
    );
}

export default GeneMapperState;