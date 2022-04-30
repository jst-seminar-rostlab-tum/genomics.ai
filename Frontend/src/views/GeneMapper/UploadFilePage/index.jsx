import { ArrowBackIcon, CheckCircleOutlinedIcon } from '@mui/icons-material/CheckCircle';
import { Button, Box, Container, Divider, Stack, Typography } from '@mui/material';
import { GeneralCard as Card } from 'components/Cards/GeneralCard';
import DatasetCard from 'components/Cards/DatasetCard';
import CustomButton from 'components/CustomButton';
import FileUpload from 'components/FileUpload';
import Input from 'components/Input/Input';
import { Modal, ModalTitle } from 'components/Modal';
import { TabGroup } from 'components/Tab';
import { useState, useEffect  } from 'react';
import { useHistory } from 'react-router-dom';
import styles from './uploadfilepage.module.css';

function UploadFilePage({ path, basePath, selectedAtlas, selectedModel, activeStep, setActiveStep, steps}) {
  const [uploadedFile, setUploadedFile] = useState();
  const [mappingName, setMappingName] = useState('');
  const [existingDatasets, setExistingDatasets] = useState([]);
  const [ongoingUploads, setOngoingUploads] = useState([]);
  const [tabsValue, setTabsValue] = useState(0);
  const [requirements, setRequirements] = useState([]);
  const [open, setOpen] = useState(false);
  const [atlasInfoOpen, setAtlasInfoOpen] = useState(false);
  const [modelInfoOpen, setModelInfoOpen] = useState(false);
  const history = useHistory();

  const [tabLabels] = useState([
    {
      label: 'Exisiting Datasets',
      additionalContent: (
        <Box sx={{ flexDirection:'column', maxHeight: '10.5em'}}>
          <DatasetCard title='dataset1' width='96%' height='3em'/>
          <DatasetCard title='dataset2' width='96%' height='3em'/>
          <DatasetCard title='dataset3' width='96%' height='3em'/>
          <DatasetCard title='dataset4' width='96%' height='3em'/>
          <DatasetCard title='dataset5' width='96%' height='3em'/>
          {/* { existingDatasets ? existingDatasets.map(data => {
              <DatasetCard title='test' width='95%' height='3em' />
            }) : 
            <Typography>No existing datasets available.</Typography> 
          } */}
        </Box>
      )
    },
    {
      label: 'Ongoing Uploads',
      additionalContent: (
        <Box sx={{ flexDirection:'column', maxHeight: '10.5em'}}>
          <DatasetCard title='upload1' width='96%' height='3em'/>
          <DatasetCard title='upload2' width='96%' height='3em'/>
          <DatasetCard title='upload3' width='96%' height='3em'/>
          <DatasetCard title='upload4' width='96%' height='3em'/>
          <DatasetCard title='upload5' width='96%' height='3em'/>
          {/* { ongoingUploads ? ongoingUploads.map((data, index) => {
              <Card title={data} category='in progress' />
            }) : 
            <Typography>There are currently no ongoing uploads.</Typography> 
          } */}
        </Box>
      )
    },
  ]);

  useEffect(() => {
    // TODO make API calls here for existing datasets and ongoing uploads
  })

  const handleOnDropChange = (file) => {
    console.log(file)
    setUploadedFile(file);
  }

  const handleSubmit = () => {
    // save mapping name
    // TODO: File validation on frontend first? (only one file could be uploaded at once, etc.)
    // TODO: POST request
    setOpen(false); // opens modal to input mapping name
    history.push(`${path}`);  // go back to GeneMapper home
  }

  return (
    <Container>
      <Stack
        direction="row"
        divider={(<Divider className={styles.divider}orientation="vertical" flexItem />)}
      >
      {/* Left side */}
      <Container>
        <Stack direction="column">
          <Typography sx={{ fontWeight: 'bold', fontSize: '20px', paddingBottom:'1em' }}>Your Choice</Typography>
          <Stack direction='row' spacing={2} sx={{ paddingBottom:'1.5em' }}>
            <Card
              width='50%'
              children={(
                <Stack direction='column'>
                  {/* TODO atlas details */}
                  <Typography sx={{ fontWeight: 'bold', fontSize: '18px'}}>{selectedAtlas ? selectedAtlas : "Atlas"}</Typography>
                  <Typography sx={{ fontSize: '12px'}}>
                    Lorem ipsum dolor sit amet, consetetur sadipscing elitr, sed diam nonumy eirmod tempor  
                  </Typography>
                  <Button size="small" disabled={!selectedAtlas} onClick={() => setAtlasInfoOpen(true)}>Learn More</Button>
                  <Modal isOpen={atlasInfoOpen} setOpen={setAtlasInfoOpen} children={(
                    <Container>
                      <ModalTitle>{selectedAtlas}</ModalTitle>
                      <Typography>{/* Atlas information here */}</Typography>
                      <Button size="large" onClick={() => setAtlasInfoOpen(false)}>Close</Button>
                    </Container>
                    )}
                  />
                </Stack>
                )}
            />
            <Card
              width='50%'
              children={(
                <Stack direction='column'>
                  {/* TODO model details */}
                  <Typography sx={{ fontWeight: 'bold', fontSize: '18px'}}>{ selectedModel ? selectedModel : "Model"}</Typography>
                  <Typography sx={{ fontSize: '12px'}}>
                  Lorem ipsum dolor sit amet, consetetur sadipscing elitr, sed diam nonumy eirmod tempor  
                  </Typography>
                  <Button size="small" disabled={!selectedModel} onClick={() => setModelInfoOpen(true)}>Learn More</Button>
                  <Modal isOpen={modelInfoOpen} setOpen={setModelInfoOpen} children={(
                    <Container>
                      <ModalTitle>{selectedModel}</ModalTitle>
                      <Typography>{/* Model information here */}</Typography>
                      <Button size="large" onClick={() => setModelInfoOpen(false)}>Close</Button>
                    </Container>
                    )}
                  />
                </Stack>
              )}
            />
          </Stack>
          <Card
            children={(
              // TODO: card styling here?
                <Stack direction='column'>
                  <Typography sx={{ fontWeight: 'bold', fontSize: '18px'}}>Consquent Requirements</Typography>
                  { requirements.length !== 0 ? 
                    requirements.map((text) => (
                      <Typography sx={{ fontSize: '12px' }}>
                        <li> {text} </li>
                      </Typography>
                    )) : 
                    <Typography variant='h7'>
                      There are no consequent requirements!
                    </Typography>
                  }
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
              <ModalTitle> Give your mapping a name </ModalTitle>
              <Divider className={styles.divider} orientation="horizontal" flexItem />
              <Input 
                placeholder="Enter name here"
                defaultValue={mappingName}
                isRequired
              />
              <Stack direction='row'>
                <Button size="large" onClick={() => setOpen(false)}>Close</Button>
                <Button size="large" onClick={handleSubmit}>Done</Button>
              </Stack>
            </Container>
            )}
          />
          <Stack className="flexContainer" direction="column">
            <Typography sx={{ fontWeight: 'bold', fontSize: '20px', paddingBottom:'1.5em' }}>Upload Datasets</Typography>
            <FileUpload 
              height='200px'
              handleFileUpload={handleOnDropChange}
            />
          </Stack>
          <Box>
            <TabGroup tabsInfo={tabLabels} value={tabsValue} setValue={setTabsValue} />
          </Box>
        </Container>
      </Stack>
      <Stack direction='row' justifyContent='space-between' sx={{ marginTop:'75px'}}>
        <CustomButton type='tertiary' onClick={() => setActiveStep(0)}>
          Back
        </CustomButton>
        <CustomButton type='primary' onClick={() => {
          setOpen(true);
        }}>
          Create Mapping
        </CustomButton>
      </Stack>
    </Container>
  );
}

export default UploadFilePage;
