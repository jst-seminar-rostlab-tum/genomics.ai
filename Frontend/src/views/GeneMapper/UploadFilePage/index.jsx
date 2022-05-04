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
import React, { useState, useEffect, useCallback } from 'react';
import { useHistory } from 'react-router-dom';
import styles from './uploadfilepage.module.css';
import ProjectMock from 'shared/services/mock/projects';
import ProjectService from 'shared/services/Project.service';
import { useSubmissionProgress } from 'shared/context/submissionProgressContext';
import { TabCard } from 'components/GeneMapper/TabCard';

function UploadFilePage({
  path, selectedAtlas, selectedModel, setActiveStep,
}) {
  const [uploadedFile, setUploadedFile] = useState();
  const [selectedDataset, setSelectedDataset] = useState();
  const [mappingName, setMappingName] = useState('');
  const [existingDatasets, setExistingDatasets] = useState([]);
  const [requirements, setRequirements] = useState([]);
  const [open, setOpen] = useState(false);
  const [atlasInfoOpen, setAtlasInfoOpen] = useState(false);
  const [modelInfoOpen, setModelInfoOpen] = useState(false);
  const [submissionProgress, setSubmissionProgress] = useSubmissionProgress();
  const history = useHistory();

  useEffect(() => {
    setRequirements(selectedModel.requirements || [
      // source: https://beta.fastgenomics.org/analyses/scarches
      'Ensure your data is in h5ad format',
      'Make sure your .X layer has raw counts (i.e. integers, so no normalization, no log-transformation)',
      'If your dataset contains multiple batches, specify these in the .obs layer under .obs["dataset"]',
    ]);
  }, [selectedModel]);

  useEffect(() => {
    ProjectMock.getDatasets().then((data) => {
      setExistingDatasets(data);
    });
  }, [existingDatasets]);

  const handleOnDropChange = (file) => {
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
    console.log(selectedDataset);
    // save mapping name
    setOpen(false); // opens modal to input mapping name
    createProject(mappingName, selectedAtlas._id, selectedModel._id, uploadedFile ? uploadedFile[0] : selectedDataset);
    history.push(`${path}`); // go back to GeneMapper home
  };

  const handleSelectDataset = (data) => {
    if (data.name === selectedDataset?.name) {
      setSelectedDataset(null);
    } else {
      setSelectedDataset(data);
    }
  };

  return (
    <Box sx={{ marginTop: '2.5em' }}>
      <Stack
        direction="row"
        divider={(<Divider className={styles.divider} orientation="vertical" flexItem />)}
      >
        {/* Left side */}
        <Box width="50%" mr="3%">
          <Stack direction="column">
            <Typography variant="h5" fontWeight="bold" pb="1em">Your Choice</Typography>
            <Stack direction="row" spacing={2} sx={{ paddingBottom: '1.5em' }}>
              <Card
                width="50%"
                children={(
                  <Stack direction="column">
                    <Typography gutterBottom variant="h6" fontWeight="bold">{selectedAtlas.name}</Typography>
                    <Typography gutterBottom variant="caption">
                      { `Species: ${selectedAtlas.species}` }
                    </Typography>
                    <Button size="small" disabled={!selectedAtlas} onClick={() => setAtlasInfoOpen(true)}>Learn More</Button>
                    <Modal
                      isOpen={atlasInfoOpen}
                      setOpen={setAtlasInfoOpen}
                      children={(
                        <Container>
                          <ModalTitle>{selectedAtlas.name}</ModalTitle>
                          <Typography variant="body1" gutterBottom>
                            {
                              /* Atlas information here */
                              Object.keys(selectedAtlas).map((key, i) => (<li key={i}>{`${key} : ${selectedAtlas[key]}`}</li>))
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
                    <Typography gutterBottom variant="h6" fontWeight="bold">{selectedModel.name}</Typography>
                    <Typography
                      gutterBottom
                      variant="caption"
                      sx={{
                        whiteSpace: 'nowrap', overflow: 'hidden', textOverflow: 'ellipsis', maxWidth: '200px',
                      }}
                    >
                      {selectedModel.description}
                    </Typography>
                    <Button size="small" disabled={!selectedModel} onClick={() => setModelInfoOpen(true)}>Learn More</Button>
                    <Modal
                      isOpen={modelInfoOpen}
                      setOpen={setModelInfoOpen}
                      children={(
                        <Container>
                          <ModalTitle>{selectedModel.name}</ModalTitle>
                          <Typography variant="body1" gutterBottom>
                            {
                              /* Model information here */
                              Object.keys(selectedModel).map((key, i) => (<li key={i}>{`${key} : ${selectedModel[key]}`}</li>))
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
            <Stack>
              <Typography variant="h5" fontWeight="bold" pb="1em">Consequent Requirements</Typography>
              <Card>
                <Box sx={{ flexDirection: 'column', minHeight: '18em' }}>
                  {requirements
                    ? requirements.map((text) => (
                      <Box key={text} sx={{ display: 'flex' }}>
                        <Typography variant="body2" sx={{ pr: 1, fontWeight: 'bold' }}>-</Typography>
                        <Typography variant="body2">
                          {text}
                        </Typography>
                      </Box>
                    ))
                    : (
                      <Typography variant="body2" gutterBottom>
                        There are no consequent requirements!
                      </Typography>
                    )}
                </Box>

              </Card>
            </Stack>
          </Stack>
        </Box>
        {/* Right side */}
        <Box width="50%" ml="3%">
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
                  isRequired
                />
                <Stack direction="row">
                  <Button size="large" onClick={() => setOpen(false)}>Close</Button>
                  <Button size="large" onClick={handleSubmit}>Done</Button>
                </Stack>
              </Container>
            )}
          />
          <Stack className="flexContainer" direction="column" pb="1em">
            <Typography variant="h5" fontWeight="bold" pb="1em">Select Datasets for Upload</Typography>
            <FileUpload
              height="250px"
              handleFileUpload={handleOnDropChange}
            />
          </Stack>
          <Stack mt="1em" maxHeight="50%">
            <Typography variant="h5" fontWeight="bold" pb="1em">Select Existing Datasets</Typography>
            { existingDatasets ? existingDatasets.map((data) => <TabCard data={data} width="95%" height="3em" handleOnClick={() => handleSelectDataset(data)} selected={selectedDataset && data._id === selectedDataset._id} />)
              : <Typography>No existing datasets available.</Typography>}
          </Stack>
        </Box>
      </Stack>
      <Stack direction="row" justifyContent="space-between" sx={{ marginTop: '75px', marginBottom: '3em' }}>
        <CustomButton type="tertiary" onClick={() => setActiveStep(0)}>
          <ArrowBackIcon sx={{ marginRight: '2px' }} />
          Back
        </CustomButton>
        <Stack direction="row" spacing={3} alignItems="center">
          <Typography variant="h6" fontWeight="bold">{uploadedFile && `Selected file: ${uploadedFile[0].name}`}</Typography>
          <CustomButton
            type="primary"
            disabled={!uploadedFile && !selectedDataset}
            onClick={() => {
              setOpen(true);
            }}
          >
            Create Mapping
            <CheckCircleOutlineIcon sx={{ marginLeft: '4px' }} />
          </CustomButton>
        </Stack>
      </Stack>
    </Box>
  );
}

export default UploadFilePage;
