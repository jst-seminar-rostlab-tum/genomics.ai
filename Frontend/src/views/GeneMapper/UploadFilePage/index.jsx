import { ArrowBackIcon, CheckCircleOutlinedIcon } from '@mui/icons-material/CheckCircle';
import { Box, Button, Container, Divider, Stack, Step, StepLabel, Stepper, Typography } from '@mui/material';
import { GeneralCard as Card } from 'components/Cards/GeneralCard';
import CustomButton from 'components/CustomButton';
import FileUpload from 'components/FileUpload';
import Input from 'components/Input/Input';
import { Modal, ModalTitle } from 'components/Modal';
import { TabGroup } from 'components/Tab';
import { useState } from 'react';
import styles from './uploadfilepage.module.css';

function UploadFilePage({ basePath, selectedAtlas, selectedModel, activeStep, setActiveStep, steps}) {
  const [uploadedFile, setUploadedFile] = useState();
  const [existingDatasets, setExistingDatasets] = useState([]);
  const [ongoingUploads, setOngoingUploads] = useState([]);
  const [tabsValue, setTabsValue] = useState(0);
  const [requirements, setRequirements] = useState([]);
  const [open, setOpen] = useState(false);
  // attribute to go back

  const [tabLabels] = useState([
    {
      label: 'Exisiting Datasets',
      path: `${basePath}/datasets`,
    },
    {
      label: 'Ongoing Uploads',
      path: `${basePath}/ongoing`,
    },
  ]);

  const handleOnDropChange = (file) => {
    console.log(file)
    setUploadedFile(file);
  }

  const handleSubmit = () => {
    // TODO: File validation on frontend first? (only one file could be uploaded at once, etc.)
    // TODO: POST request
    setOpen(false); // opens modal to input mapping name
  }

  return (
    <Container>
      <Stack
        direction="row"
        divider={(
          <Divider 
            className={styles.divider}
            orientation="vertical"
            flexItem
          />
        )}
      >
      {/* Left side */}
      <Container>
        <Stack direction="column">
          <Typography sx={{ fontWeight: 'bold', fontSize: '20px', paddingBottom:'1em' }}>Your Choice</Typography>
          <Stack direction='row' spacing={2} sx={{ paddingBottom:'1.5em' }}>
            <Card
              width='50%'
              component={(
                <Stack direction='column'>
                  <Typography sx={{ fontWeight: 'bold', fontSize: '18px'}}>{selectedAtlas}</Typography>
                  <Typography sx={{ fontSize: '12px'}}>
                    Lorem ipsum dolor sit amet, consetetur sadipscing elitr, sed diam nonumy eirmod tempor  
                  </Typography>
                  <Button size="small">Learn More</Button>
                </Stack>
                )}
            />
            <Card
              width='50%'
              component={(
                <Stack direction='column'>
                  <Typography sx={{ fontWeight: 'bold', fontSize: '18px'}}>Model</Typography>
                  <Typography sx={{ fontSize: '12px'}}>
                  Lorem ipsum dolor sit amet, consetetur sadipscing elitr, sed diam nonumy eirmod tempor  
                  </Typography>
                  <Button size="small">Learn More</Button>
                </Stack>
              )}
            />
          </Stack>
          <Card
            component={(
              // TODO: card styling here?
                <Stack direction='column'>
                  <Typography sx={{ fontWeight: 'bold', fontSize: '18px'}}>Consquent Requirements</Typography>
                  { requirements.length !== 0 ? 
                    requirements.map((text, index) => (
                      <Typography sx={{ fontSize: '12px' }}> {++index + ' ' + text } </Typography>
                    )) : 
                    <Typography sx={{ fontSize: '14px' }}>There are no consequent requirements!</Typography>}
                  <Typography sx={{ fontSize: '12px'}}>
                    Lorem ipsum dolor sit amet, consetetur sadipscing elitr, sed diam nonumy eirmod tempLorem ipsum dolor sit amet, consetetur sadipscing elitr, sed diam
                    Lorem ipsum dolor sit amet, consetetur
                    Lorem ipsum dolor sit amet, consetetur sadipscing elitr, sed diam nonumy eirmod tempor invidunt ut labore et dolore magna aliquyam
                    Lorem ipsum dolor sit amet, consetetur sadipscing elitr, sed diam nonumy eirmod tempor invidunt ut labore et dolore magna aliquyam erat, sed diam voluptua. At vero eos et accusam et
                    Lorem ipsum dolor sit amet, conseteturor  
                  </Typography>
                </Stack>
            )}
          />
        </Stack>
        </Container>
        {/* Right side */}
        <Container>
          <Modal isOpen={open} setOpen={setOpen} children={(
            <Container>
              <ModalTitle 
                  children={(<Divider className={styles.divider}orientation="vertical" flexItem
                />)}
              >Give your mapping a name:</ModalTitle>
              <Input placeholder="Enter name here" />
              <Stack direction='row'>
                <Button size="large" onClick={() => setOpen(false)}>Close</Button>
                <Button size="large" onClick={handleSubmit}>Done</Button>
              </Stack>
            </Container>
            )}
          />
          <Stack className="flexContainer" direction="column">
            <Typography sx={{ fontWeight: 'bold', fontSize: '20px', paddingBottom:'1em' }}>Upload Datasets</Typography>
            <FileUpload 
              height='200px'
              handleFileUpload={handleOnDropChange}
            />
          </Stack>
          <TabGroup
            darkBackground={false}
            tabsInfo={tabLabels}
            value={tabsValue}
            setValue={setTabsValue}
          >
          </TabGroup>
        </Container>
      </Stack>
      <Stack direction='row' justifyContent='space-between' sx={{ marginTop:'75px'}}>
        <CustomButton type='tertiary' onClick={() => setActiveStep(0)}>
          Back
        </CustomButton>
        <CustomButton type='primary' onClick={() => setOpen(true)} children={(
          <Stack direction='row'>
            <CheckCircleOutlinedIcon />
          </Stack> )}>
          Create Mapping
        </CustomButton>
        {/* <Button variant='outlined' startIcon={<ArrowBackIcon />}>Back</Button>
        <Button variant='contained' endIcon={<CheckCircleIcon />}>Create Mapping</Button> */}
      </Stack>
    </Container>
  );
}

export default UploadFilePage;
