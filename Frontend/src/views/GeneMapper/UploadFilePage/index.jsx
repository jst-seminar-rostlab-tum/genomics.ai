/* eslint-disable react/no-children-prop */
import ArrowBackIcon from '@mui/icons-material/ArrowBack';
import CheckCircleOutlineIcon from '@mui/icons-material/CheckCircleOutlineOutlined';
import {
  Alert, Button, Box, Container, Divider, Stack, Tooltip, Typography, TextField,
} from '@mui/material';
import { GeneralCard as Card } from 'components/Cards/GeneralCard';
import CustomButton from 'components/CustomButton';
import FileUpload from 'components/FileUpload';
import { Modal, ModalTitle } from 'components/Modal';
import React, { useState, useEffect, useCallback } from 'react';
import { useHistory, useLocation } from 'react-router-dom';
import styles from './uploadfilepage.module.css';
// import ProjectMock from 'shared/services/mock/projects';
import ProjectService from 'shared/services/Project.service';
import { initSubmissionProgress, useSubmissionProgress } from 'shared/context/submissionProgressContext';
// import { TabCard } from 'components/GeneMapper/TabCard'; 
import { LearnMoreAtlasComponent } from 'views/Explore/LearnMoreAtlas';
import { LearnMoreModelComponent } from 'views/Explore/LearnMoreModel';
import { uploadMultipartFile } from 'shared/services/UploadLogic';

function UploadFilePage({
  path, selectedAtlas, selectedModel, setActiveStep,
}) {
  const [uploadedFile, setUploadedFile] = useState();
  const [selectedDataset, setSelectedDataset] = useState();
  const [mappingName, setMappingName] = useState('');
  const [existingDatasets, setExistingDatasets] = useState();
  const [requirements, setRequirements] = useState([]);
  const [open, setOpen] = useState(false);
  const [atlasInfoOpen, setAtlasInfoOpen] = useState(false);
  const [modelInfoOpen, setModelInfoOpen] = useState(false);
  const [submissionProgress, setSubmissionProgress] = useSubmissionProgress();
  const history = useHistory();
  const { pathname } = useLocation();

  useEffect(() => {
    setRequirements(selectedModel.requirements || [
      // source: https://beta.fastgenomics.org/analyses/scarches
      'Ensure your data is in h5ad format',
      'Make sure your .X layer has raw counts (i.e. integers, so no normalization, no log-transformation)',
      'If your dataset contains multiple batches, specify these in the .obs layer under .obs["dataset"]',
    ]);
  }, [selectedModel]);

  // Temporarily commented out as the endpoint is not implemented yet
  // useEffect(() => {
  //   ProjectMock.getDatasets().then((data) => {
  //     setExistingDatasets(data);
  //   });
  // }, [existingDatasets]);

  const handleOnDropChange = (file) => {
    setUploadedFile(file);
  };
  const createProject = useCallback((projectName, atlasId, modelId, file) => {
    ProjectService.createProject(
      projectName,
      atlasId,
      modelId,
      file.name,
    ).then((project) => {
      uploadMultipartFile(
        project.uploadId,
        file,
        initSubmissionProgress(project.uploadId),
        (update) => {
          setSubmissionProgress((prev) => ({
            ...prev,
            [project._id]: update(prev[project._id] ?? initSubmissionProgress(project.uploadId)),
          }));
        },
      );
      history.push(pathname.split('/').includes('explore') ? `${path}/genemapper` : path); // go back to GeneMapper home
    });
  }, [submissionProgress]);

  const handleSubmit = (e) => {
    e?.preventDefault();
    console.log(selectedDataset);
    // save mapping name
    setOpen(false); // opens modal to input mapping name
    createProject(mappingName, selectedAtlas._id, selectedModel._id,
      uploadedFile ? uploadedFile[0] : selectedDataset);
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
                    <Typography variant="caption" fontWeight="bold">Atlas</Typography>
                    <Typography gutterBottom variant="h6" fontWeight="bold">{selectedAtlas.name}</Typography>
                    <Typography
                      gutterBottom
                      variant="caption"
                      sx={{
                        whiteSpace: 'nowrap', overflow: 'hidden', textOverflow: 'ellipsis', maxWidth: '200px',
                      }}
                    >
                      {`Modalities:  ${selectedAtlas.modalities}`}
                    </Typography>
                    <Typography gutterBottom variant="caption">{`Cells in Reference:  ${selectedAtlas.numberOfCells}`}</Typography>
                    <Typography gutterBottom variant="caption">{`Species: ${selectedAtlas.species}`}</Typography>
                    <Button
                      size="small"
                      variant="outlined"
                      onClick={() => setAtlasInfoOpen(true)}
                      sx={{
                        borderRadius: 100, width: '50%', ml: '50%', mt: '1em', textTransform: 'none',
                      }}
                    >
                      Learn more
                    </Button>
                    <Modal
                      isOpen={atlasInfoOpen}
                      setOpen={setAtlasInfoOpen}
                      children={(
                        <Container>
                          <LearnMoreAtlasComponent id={selectedAtlas._id} onClick={() => history.push(`/explore/atlases/${selectedAtlas._id}/visualization`)} />
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
                    <Typography variant="caption" fontWeight="bold">Model</Typography>
                    <Typography gutterBottom variant="h6" fontWeight="bold">{selectedModel.name}</Typography>
                    <Typography
                      gutterBottom
                      variant="caption"
                      sx={{
                        display: 'block',
                        textOverflow: 'ellipsis',
                        wordWrap: 'break-word',
                        overflow: 'hidden',
                        maxHeight: '5.7em',
                        lineHeight: '1.9em',
                      }}
                    >
                      {selectedModel.description}
                    </Typography>
                    <Button
                      size="small"
                      variant="outlined"
                      onClick={() => setModelInfoOpen(true)}
                      sx={{
                        borderRadius: 100, width: '50%', ml: '50%', mt: '1em', textTransform: 'none',
                      }}
                    >
                      Learn more
                    </Button>
                    <Modal
                      isOpen={modelInfoOpen}
                      setOpen={setModelInfoOpen}
                      children={(
                        <Container>
                          <LearnMoreModelComponent id={selectedModel._id} />
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
                    ? requirements.map((text, index) => (
                      <Box key={text} sx={{ display: 'flex' }}>
                        <Typography variant="body2" sx={{ pr: 1, fontWeight: 'bold' }}>{++index}</Typography>
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
                <form onSubmit={handleSubmit}>
                  <TextField
                    variant="standard"
                    placeholder="Enter name here"
                    fullWidth
                    onChange={(e) => setMappingName(e.target.value)}
                    required
                    label="Mapping name"
                  />
                  <Stack direction="row" justifyContent="space-between" mt="1.5em">
                    <CustomButton
                      type="tertiary"
                      onClick={() => {
                        setOpen(false);
                      }}
                    >
                      Close
                    </CustomButton>
                    <CustomButton
                      type="primary"
                      onClick={handleSubmit}
                      disabled={!mappingName}
                    >
                      Submit
                    </CustomButton>
                  </Stack>
                </form>
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
            { 
              existingDatasets ? 
                existingDatasets.map((data) => 
                  <TabCard 
                    data={data}
                    width="95%" 
                    height="3em"
                    handleOnClick={
                      () => handleSelectDataset(data)} selected={selectedDataset && data._id === selectedDataset._id
                    }
                  />)
              : <Alert severity="info"> No existing datasets available. </Alert> 
            }
          </Stack>
        </Box>
      </Stack>
      <Stack direction="row" justifyContent="space-between" sx={{ marginTop: '75px', marginBottom: '3em' }}>
        <CustomButton type="tertiary" onClick={() => setActiveStep(0)}>
          <ArrowBackIcon sx={{ marginRight: '2px' }} />
          Back
        </CustomButton>
        <Stack direction="row" spacing={3} alignItems="center">
          <Typography variant="h6" fontWeight="bold">
            { uploadedFile && `Selected file: ${uploadedFile[0].name}` }
          </Typography>
          <Tooltip title={(!uploadedFile && !selectedDataset) ? "You haven't selected or uploaded a dataset!" : ''} placement="top">
            <span>
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
            </span>
          </Tooltip>
        </Stack>
      </Stack>
    </Box>
  );
}

export default UploadFilePage;
