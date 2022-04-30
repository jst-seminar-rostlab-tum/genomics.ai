import { ArrowBackIcon, CheckCircleOutlinedIcon } from '@mui/icons-material/CheckCircle';
import {
  Button, Box, Container, Divider, Stack, Typography,
} from '@mui/material';
import { GeneralCard as Card } from 'components/Cards/GeneralCard';
import DatasetCard from 'components/Cards/DatasetCard';
import CustomButton from 'components/CustomButton';
import FileUpload from 'components/FileUpload';
import Input from 'components/Input/Input';
import { Modal, ModalTitle } from 'components/Modal';
import { TabGroup } from 'components/Tab';
import React, { useState, useEffect, useCallback } from 'react';
import { useHistory } from 'react-router-dom';
import styles from './uploadfilepage.module.css';
import ProjectMock from 'shared/services/mock/projects';
import ProjectService from 'shared/services/Project.service';
import { useSubmissionProgress } from 'shared/context/submissionProgressContext';

function UploadFilePage({
  path, selectedAtlas, selectedModel, setActiveStep
}) {
  const [uploadedFile, setUploadedFile] = useState();
  const [mappingName, setMappingName] = useState('');
  const [existingDatasets, setExistingDatasets] = useState([]);
  const [ongoingUploads, setOngoingUploads] = useState([]);
  const [tabsValue, setTabsValue] = useState(0);
  const [requirements, setRequirements] = useState([]);
  const [open, setOpen] = useState(false);
  const [atlasInfoOpen, setAtlasInfoOpen] = useState(false);
  const [modelInfoOpen, setModelInfoOpen] = useState(false);
  const [submissionProgress, setSubmissionProgress] = useSubmissionProgress();
  const history = useHistory();
  const datasets = [
    {
      _id: 1,
      name: 'dataset1',
      category: 'category1'
    },
    {
      _id: 2,
      name: 'dataset2',
      category: 'category2'
    },
    {
      _id: 3,
      name: 'dataset3',
      category: 'category3'
    },{
      _id: 4,
      name: 'dataset4',
      category: 'category4'
    },
  ]
  const [tabLabels] = useState([
    {
      label: 'Exisiting Datasets',
      additionalContent: (
        <Box sx={{ flexDirection: 'column', maxHeight: '10.5em' }}>

          {/* <DatasetCard title="dataset1" width="96%" height="3em" /> */}
          { existingDatasets.map(data => {
              console.log(data)
              return <DatasetCard title={data.name} category={data.category} width='95%' height='3em' />
            })
          }
          {/* { existingDatasets ? existingDatasets.map(data => {
              console.log(data)
              return <DatasetCard title={data.name} category={data.category} width='95%' height='3em' />
            }) :
            <Typography>No existing datasets available.</Typography>
          } */}
        </Box>
      ),
    },
    // {
    //   label: 'Ongoing Uploads',
    //   additionalContent: (
    //     <Box sx={{ flexDirection: 'column', maxHeight: '10.5em' }}>
    //       <DatasetCard title="upload1" width="96%" height="3em" />
    //       <DatasetCard title="upload2" width="96%" height="3em" />
    //       <DatasetCard title="upload3" width="96%" height="3em" />
    //       <DatasetCard title="upload4" width="96%" height="3em" />
    //       <DatasetCard title="upload5" width="96%" height="3em" />
    //       {/* { ongoingUploads ? ongoingUploads.map((data, index) => {
    //           <Card title={data} category='in progress' />
    //         }) :
    //         <Typography>There are currently no ongoing uploads.</Typography>
    //       } */}
    //     </Box>
    //   ),
    // },
  ]);

  useEffect(() => {
    setRequirements(selectedModel.requirements);
  }, []);

  useEffect(() => {
    ProjectMock.getDatasets().then((data) => setExistingDatasets(data));
  }, [existingDatasets]);

  const handleOnDropChange = (file) => {
    console.log(file);
    setUploadedFile(file);
  };

  const createProject = useCallback((projectName, atlasId, modelId, file) => {
    ProjectService.startOrContinueProjectUpload(
      file,
      submissionProgress,
      setSubmissionProgress,
      {
        projectName,
        atlasId,
        modelId,
        fileName: file.name,
      },
    );
  }, [submissionProgress]);

  const handleSubmit = () => {
    // save mapping name
    setOpen(false); // opens modal to input mapping name
    createProject(mappingName, '111122223333444455556666', '011122223333444455556666', uploadedFile[0]);
    history.push(`${path}`); // go back to GeneMapper home
  };

  return (
    <Container>
      <Stack
        direction="row"
        divider={(<Divider className={styles.divider} orientation="vertical" flexItem />)}
      >
        {/* Left side */}
        <Container>
          <Stack direction="column">
            <Typography variant='h5' sx={{ fontWeight: 'bold', paddingBottom: '1em' }}>Your Choice</Typography>
            <Stack direction="row" spacing={2} sx={{ paddingBottom: '1.5em' }}>
              <Card
                width="50%"
                children={(
                  <Stack direction="column">
                    {/* TODO atlas details */}
                    <Typography variant='h6' sx={{ fontWeight: 'bold' }}>{selectedAtlas.name}</Typography>
                    <Typography variant='caption'>
                      { 'Species: ' + selectedAtlas.species }
                    </Typography>
                    <Button size="small" disabled={!selectedAtlas} onClick={() => setAtlasInfoOpen(true)}>Learn More</Button>
                    <Modal
                      isOpen={atlasInfoOpen}
                      setOpen={setAtlasInfoOpen}
                      children={(
                        <Container>
                          <ModalTitle>{selectedAtlas.name}</ModalTitle>
                          <Typography variant='body1' gutterBottom>{
                              /* Atlas information here */
                              Object.entries(selectedAtlas).forEach(([key, value]) => <li>{key + ': ' + value}</li> )
                            }
                          </Typography>
                          <Button size="large" onClick={() => setAtlasInfoOpen(false)}>Close</Button>
                        </Container>
                    )}
                    />
                  </Stack>
                )}
              />
              <Card
                width="50%"
                children={(
                  <Stack direction="column">
                    {/* TODO model details */}
                    <Typography variant='h6' gutterBottom sx={{ fontWeight: 'bold' }}>{selectedModel.name}</Typography>
                    <Typography variant='caption' gutterBottom>
                      {selectedModel.description}
                    </Typography>
                    <Button size="small" disabled={!selectedModel} onClick={() => setModelInfoOpen(true)}>Learn More</Button>
                    <Modal
                      isOpen={modelInfoOpen}
                      setOpen={setModelInfoOpen}
                      children={(
                        <Container>
                          <ModalTitle>{selectedModel.name}</ModalTitle>
                          <Typography variant='body1' gutterBottom>{
                              /* Model information here */
                              Object.entries(selectedModel).forEach(([key, value]) => <li>{key + ': ' + value}</li> )
                            }
                          </Typography>
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
                <Box sx={{ flexDirection:'column', backgroundColor: 'e9ecef'}}>
                  <Typography variant='h6' sx={{ fontWeight: 'bold' }}>Consquent Requirements</Typography>
                  { requirements
                    ? requirements.map((text) => (
                      <Typography variant='body2' gutterBottom>
                        <li>{text}</li>
                      </Typography>
                    ))
                    : (
                      <Typography variant="body2" gutterBottom sx={{ fontWeight:'bold' }}>
                        There are no consequent requirements!
                      </Typography>
                    )}
                </Box>
            )}
            />
          </Stack>
        </Container>
        {/* Right side */}
        <Container>
          <Modal
            isOpen={open}
            setOpen={setOpen}
            children={(
              <Container>
                <ModalTitle> Give your mapping a name </ModalTitle>
                <Divider className={styles.divider} orientation="horizontal" flexItem />
                <Input
                  placeholder="Enter name here"
                  defaultValue={mappingName}
                  onChange={(e) => setMappingName(e.target.value)}
                  isRequired
                />
                <Stack direction="row">
                  <Button size="large" onClick={() => setOpen(false)}>Close</Button>
                  <Button size="large" onClick={handleSubmit}>Done</Button>
                </Stack>
              </Container>
            )}
          />
          <Stack className="flexContainer" direction="column">
            <Typography variant='h5' sx={{ fontWeight: 'bold', paddingBottom: '1em' }}>Upload Datasets</Typography>
            <FileUpload
              height="200px"
              handleFileUpload={handleOnDropChange}
            />
          </Stack>
          <Box>
            <TabGroup tabsInfo={tabLabels} value={tabsValue} setValue={setTabsValue} />
          </Box>
        </Container>
      </Stack>
      <Stack direction="row" justifyContent="space-between" sx={{ marginTop: '75px' }}>
        <CustomButton type="tertiary" onClick={() => setActiveStep(0)}>
          Back
        </CustomButton>
        <CustomButton
          type="primary"
          onClick={() => {
            setOpen(true);
          }}
        >
          Create Mapping
        </CustomButton>
      </Stack>
    </Container>
  );
}

export default UploadFilePage;
