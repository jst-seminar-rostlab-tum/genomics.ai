/* eslint-disable react/no-children-prop */
import ArrowBackIcon from '@mui/icons-material/ArrowBack';
import CheckCircleOutlineIcon from '@mui/icons-material/CheckCircleOutlineOutlined';
import {
  Alert, Box, Button, Container, Divider, Stack, TextField, Tooltip, Typography,
} from '@mui/material';
import { GeneralCard as Card } from 'components/Cards/GeneralCard';
import CustomButton from 'components/CustomButton';
import FileUpload from 'components/FileUpload';
import { Modal, ModalTitle } from 'components/Modal';
import { useCallback, useEffect, useState } from 'react';
import { DemoDatasetCard } from 'components/GeneMapper/TabCard';
import { useHistory } from 'react-router-dom';
import { initSubmissionProgress, useSubmissionProgress } from 'shared/context/submissionProgressContext';
import ProjectService from 'shared/services/Project.service';
import { uploadMultipartFile } from 'shared/services/UploadLogic';
import { LearnMoreAtlasComponent } from 'views/Explore/LearnMoreAtlas';
import { LearnMoreModelComponent } from 'views/Explore/LearnMoreModel';
import styles from './uploadfilepage.module.css';

function UploadFilePage({
  path, selectedAtlas, selectedModel, setActiveStep,
}) {
  const demoDatasets = [
    {
      demoId: 12,
      title: 'Pancreas + SCVI',
      atlas: 'Pancreas',
      model: 'SCVI',
      url: 'link to the dataset to fetch',
    },
    {
      demoId: 123,
      title: 'pancreas + SCANVI',
      atlas: 'Pancreas',
      model: 'SCANVI',
      url: 'link to the dataset to fetch',
    },
    {
      demoId: 1231,
      title: 'fetal + SCVI',
      atlas: 'fetal immune atlas',
      model: 'SCANVI',
      url: 'link to the dataset to fetch',
    },
    {
      demoId: 123131,
      title: 'pbmc + totalvi',
      atlas: 'pbmc',
      model: 'totalvi',
      url: 'link to the dataset to fetch',
    },
    { // Leave this object in since it is causing bugs. Find out why
      demoId: 131131,
      title: 'Fetal Immune + TotalVI',
      atlas: 'Fetal Immune',
      model: 'TOTALVI',
      url: 'link to the dataset to fetch',
    },
  ];

  const [uploadedFile, setUploadedFile] = useState();
  const [selectedDataset, setSelectedDataset] = useState();
  const [mappingName, setMappingName] = useState('');
  const [existingDatasets, setExistingDatasets] = useState(demoDatasets);
  const [requirements, setRequirements] = useState([]);
  const [open, setOpen] = useState(false);
  const [atlasInfoOpen, setAtlasInfoOpen] = useState(false);
  const [modelInfoOpen, setModelInfoOpen] = useState(false);
  const [submissionProgress, setSubmissionProgress] = useSubmissionProgress();
  const [showWarning, setShowWarning] = useState(false);
  const [showFileWarning, setShowFileWarning] = useState(false);
  const [showAcceptedFile, setShowAcceptedFile] = useState(false);
  const history = useHistory();

  useEffect(() => {
    setRequirements(selectedModel.requirements);
  }, [selectedModel]);

  const handleOnDropChange = (file) => {
    setUploadedFile(file);
  };

  // custom file extension validator
  const validateUploadFile = (file) => {
    if (file.name.split('.').pop() !== 'h5ad') {
      setShowFileWarning(true);
      // error object returned in case of rejection
      return {
        code: 'Wrong file format',
        message: 'File must be in h5ad format',
      };
    }
    setShowFileWarning(false);

    return null; // file accepted
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
      history.push(path); // go back to GeneMapper home
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
      {showWarning
        && (
        <Alert severity="error" xs={{ marginTop: '100px' }}>
          Select or upload a dataset before continuing
        </Alert>
        )}
      <Stack
        direction="row"
        divider={(<Divider className={styles.divider} orientation="vertical" flexItem />)}
        sx={showWarning ? { marginTop: '1em' } : {}}
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
                <Box sx={{ flexDirection: 'column', minHeight: '6em' }}>
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
              validator={validateUploadFile}
              rejectionHandler={() => setUploadedFile()}
            />
            {
              showFileWarning
              && (
              <Alert severity="error" sx={{ marginTop: '1em' }}>
                File must be in h5ad format.
              </Alert>
              )
            }
            {
              uploadedFile && uploadedFile[0]
              && (
              <Alert severity="success" sx={{ marginTop: '1em' }}>
                Selected file:
                {' '}
                {uploadedFile[0].name}
              </Alert>
              )
            }
          </Stack>
          <Stack maxHeight="50%">
            <Typography variant="h5" fontWeight="bold" mb="0.5em">Select Existing Datasets</Typography>
            {
              // demoDatasets.length !==0 
              //   ? 
                [1,2,3,4].map((data) => (
                  <DemoDatasetCard
                    data={data}
                    width="95%"
                    height="3em"
                    handleOnClick={
                      () => handleSelectDataset(data)
}
                    selected={selectedDataset && data._id === selectedDataset._id}
                  />
                ))
                // : <Alert severity="info"> No existing datasets available. </Alert>
            }
          </Stack>
        </Box>
      </Stack>
      <Stack direction="row" justifyContent="space-between" sx={{ marginTop: '2em', marginBottom: '3em' }}>
        <CustomButton type="tertiary" onClick={() => setActiveStep(0)}>
          <ArrowBackIcon sx={{ marginRight: '2px' }} />
          Back
        </CustomButton>
        <Stack direction="row" spacing={3} alignItems="center">
          <Typography variant="h6" fontWeight="bold" />
          <Tooltip title={(!uploadedFile && !selectedDataset) ? "You haven't selected or uploaded a dataset!" : ''} placement="top">
            <Box onClick={!uploadedFile && !selectedDataset ? () => setShowWarning(true) : () => {}}>
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
            </Box>
          </Tooltip>
        </Stack>
      </Stack>
    </Box>
  );
}

export default UploadFilePage;
