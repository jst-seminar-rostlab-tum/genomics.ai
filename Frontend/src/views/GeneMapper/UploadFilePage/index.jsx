import ArrowBackIcon from '@mui/icons-material/ArrowBack';
import CheckCircleOutlineIcon from '@mui/icons-material/CheckCircleOutlineOutlined';
import {
  Button, Box, Container, Divider, Stack, Typography,
} from '@mui/material';
import { GeneralCard as Card } from 'components/Cards/GeneralCard';
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
import { TabCard } from 'components/GeneMapper/TabCard';
import { colors } from 'shared/theme/colors';

function UploadFilePage({
  path, selectedAtlas, selectedModel, setActiveStep
}) {
  const [uploadedFile, setUploadedFile] = useState();
  const [mappingName, setMappingName] = useState('');
  const [existingDatasets, setExistingDatasets] = useState([]);
  // const [ongoingUploads, setOngoingUploads] = useState([]);
  const [tabsValue, setTabsValue] = useState(0);
  const [requirements, setRequirements] = useState([]);
  const [open, setOpen] = useState(false);
  const [atlasInfoOpen, setAtlasInfoOpen] = useState(false);
  const [modelInfoOpen, setModelInfoOpen] = useState(false);
  const [openUploadConfirmation, setOpenUploadConfirmation] = useState(false);
  const [submissionProgress, setSubmissionProgress] = useSubmissionProgress();
  const history = useHistory();

  const [tabLabels, setTabLabels] = useState([
    {
      label: 'Exisiting Datasets',
      additionalContent: (
        <Typography>Loading datasets...</Typography>
      ),
    },
  ]);

  useEffect(() => {
    setRequirements(selectedModel.requirements);
  }, []);

  useEffect(() => {
    ProjectMock.getDatasets().then((data) => {
      setExistingDatasets(data);
      setTabLabels(
        [
          {
            label: 'Exisiting Datasets',
            additionalContent: (
              <Box sx={{ flexDirection: 'column', maxHeight: '50%' }}>
                { existingDatasets ? existingDatasets.map(data => {
                    return <TabCard fileName={data.name} status={data.status} width='95%' height='3em' />
                  }) :
                  <Typography>No existing datasets available.</Typography>
                }
              </Box>
            ),
          },
        ]
      )
    });
  }, [existingDatasets]);

  const handleOnDropChange = (file) => {
    // console.log(file);
    setUploadedFile(file);
    setOpenUploadConfirmation(true);
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
    <Box sx={{ marginTop: '2.5em' }}>
      <Stack
        direction="row"
        divider={(<Divider className={styles.divider} orientation="vertical" flexItem />)}
        justifyContent="space-between"
      >
        {/* Left side */}
        <Container>
          <Stack direction="column">
            <Typography variant='h5' fontWeight='bold' pb='1em'>Your Choice</Typography>
            <Stack direction="row" spacing={2} sx={{ paddingBottom: '1.5em' }}>
              <Card
                width="50%"
                children={(
                  <Stack direction="column">
                    <Typography gutterBottom variant='h6' fontWeight='bold'>{selectedAtlas.name}</Typography>
                    <Typography gutterBottom variant='caption'>
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
                              Object.keys(selectedAtlas).map((key, i) => {
                                return (<li key={i}>{key + ' : ' + selectedAtlas[key]}</li>)
                              })
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
                    <Typography gutterBottom variant='h6' fontWeight='bold'>{selectedModel.name}</Typography>
                    <Typography gutterBottom variant='caption' sx={{ whiteSpace: 'nowrap', overflow: 'hidden', textOverflow: 'ellipsis', maxWidth: '200px'}}>
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
                              Object.keys(selectedModel).map((key, i) => {
                                return (<li key={i}>{key + ' : ' + selectedModel[key]}</li>)
                              })
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
              bg={colors.neutral[100]}
              children={(
                <Box sx={{ flexDirection:'column', minHeight: '19.5em'}}>
                  <Typography variant='h6' fontWeight='bold'>Consquent Requirements</Typography>
                  { requirements
                    ? requirements.map((text) => (
                      <Typography variant='body2' gutterBottom>
                        <li>{text}</li>
                      </Typography>
                    ))
                    : (
                      <Typography variant="body2" gutterBottom fontWeight='bold'>
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
                <ModalTitle>Give your mapping a name </ModalTitle>
                <Divider className={styles.divider} orientation="horizontal" flexItem />
                <Input
                  placeholder="Enter name here"
                  defaultValue={mappingName}
                  onChangeEvent={setMappingName}
                  disabledHandler
                  isRequired
                />
                <Stack direction="row">
                  <Button size="large" onClick={() => setOpen(false)}>Close</Button>
                  <Button size="large" onClick={handleSubmit}>Done</Button>
                </Stack>
              </Container>
            )}
          />
          {/* Modal to confirm that file is being uploaded */}
          <Modal
            isOpen={openUploadConfirmation}
            setOpen={setOpenUploadConfirmation}
            children={(
              <Container>
                <ModalTitle>File Upload Confirmation</ModalTitle>
                <Typography variant='body1' gutterBottom>
                  { uploadedFile && 'The file \"' +  uploadedFile[0].name + '\" is being uploaded!'}
                </Typography>
                <Button size="large" onClick={() => setOpenUploadConfirmation(false)}>Close</Button>
              </Container>
            )}
          />
          <Stack className="flexContainer" direction="column">
            <Typography variant='h5' fontWeight='bold' pb='1em'>Upload Datasets</Typography>
            <FileUpload
              height="250px"
              handleFileUpload={handleOnDropChange}
            />
          </Stack>
          <Box>
            <TabGroup tabsInfo={tabLabels} value={tabsValue} setValue={setTabsValue} />
          </Box>
        </Container>
      </Stack>
      <Stack direction="row" justifyContent="space-between" sx={{ marginTop: '75px', marginBottom: '3em'}}>
        <CustomButton type="tertiary" onClick={() => setActiveStep(0)}>
          <ArrowBackIcon sx={{marginRight: '2px'}} />
          Back
        </CustomButton>
        <CustomButton
          type="primary"
          onClick={() => {
            setOpen(true);
          }}
        >
          Create Mapping
          <CheckCircleOutlineIcon sx={{marginLeft: '4px'}} />
        </CustomButton>
      </Stack>
    </Box>
  );
}

export default UploadFilePage;
