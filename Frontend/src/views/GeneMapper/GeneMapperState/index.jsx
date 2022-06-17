import {
  Box, Container, Step, StepButton, Stepper,
} from '@mui/material';
import React, { useState, useEffect } from 'react';
import AtlasModelChoice from '../AtlasModelChoice/AtlasModelChoice';
import UploadFilePage from '../UploadFilePage';
import { useLocation, useHistory } from 'react-router-dom';
import ModelService from 'shared/services/Model.service';
import AtlasService from 'shared/services/Atlas.service';

function GeneMapperState({ path }) {
  const { search } = useLocation();
  const history = useHistory();
  const searchParams = new URLSearchParams(search);
  const urlParams = new URLSearchParams(searchParams);
  const filterParams = Object.fromEntries(urlParams);
  // function to fetch the demo datasets, currently an object but
  const demoDatasetArr = [
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
      demoId: 12312,
      title: 'fetal + SCVI 2',
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
  ];

  const atlasId = filterParams.atlas;
  const modelId = filterParams.model;
  const [selectedAtlas, setSelectedAtlas] = useState('');
  const [selectedModel, setSelectedModel] = useState('');
  const [activeStep, setActiveStep] = useState(0);
  const [selectedDemoDataset, setSelectedDemoDataset] = useState(null);
  const [demoDatasetIsSelected, setDemoDatasetIsSelected] = useState(false);
  const [demoDatasets, setDemoDatasets] = useState(demoDatasetArr);
  const steps = ['Pick Atlas and Model', 'Choose File and Project details'];
  const [atlases, setAtlases] = useState(null);
  const [models, setModels] = useState(null);

  const handleAtlasSelection = (newAtlas) => {
    setSelectedAtlas(newAtlas);
    setSelectedModel('');
  };

  // todo: function that fetches the demo datasets

  // function to update the state in the URL
  const updateQueryParams = (param, value) => {
    const params = new URLSearchParams(history.location.search);
    if (value) {
      params.set(param, value);
    } else {
      params.delete(param);
    }

    history.replace({
      pathname: history.location.pathname,
      search: params.toString(),
    });
  };

  const handleStep = (step) => {
    setActiveStep(step);
    if (step === 1 && selectedAtlas && selectedModel) {
      updateQueryParams('atlas', selectedAtlas._id);
      updateQueryParams('model', selectedModel._id);
    }
  };

  useEffect(() => {
    AtlasService.getAtlases().then((a) => {
      a.map((a) => {
        // adjust numberOfCells
        let { numberOfCells } = a;
        let dimension = '';
        if (numberOfCells > 1000000000) {
          numberOfCells = Math.round(numberOfCells / 1000000000);
          dimension = 'B';
        } else if (numberOfCells > 1000000) {
          numberOfCells = Math.round(numberOfCells / 1000000);
          dimension = 'M';
        } else if (numberOfCells > 1000) {
          numberOfCells = Math.round(numberOfCells / 1000);
          dimension = 'K';
        }
        a.numberOfCells = numberOfCells + dimension;

        // adjust modalities
        if (!(typeof a.modalities === 'string')) {
          // modalities ist array of strings
          if (a.modalities.length == 0) {
            a.modalities = 'None';
          } else if (a.modalities.length == 1) {
            a.modalities = a.modalities[0];
          } else {
            a.modalities = `${a.modalities[0]}, ...`;
          }
        }
        if (a._id === atlasId) {
          setSelectedAtlas(a);
        }
        return a;
      });
      setAtlases(a);
    });

    ModelService.getModels().then((m) => {
      m.map((model) => {
                model.requirements = [
                    'Ensure your data is in h5ad format',
                    'Batch/study information is mandatory and should be labeled as “batch”'
                ];
        if (model.name === 'scVI') {
          model.requirements.push('Cell type information should be labeled as “cell_type”');
          model.requirements.push('For unlabeled cells, the value for “cell_type” should be “Unknown”');
        }
        if (model.name === 'scANVI') {
          model.requirements.push('Cell type information should be labeled as “cell_type”');
        }
        if (model._id === modelId) {
          setSelectedModel(model);
        }
      });
      setModels(m);
    });
  }, []);

  useEffect(() => {
    if (atlasId && selectedAtlas && modelId && selectedModel) {
      handleStep(1);
    } else if (atlases && models) {
      history.push({
        pathname: history.location.pathname,
        search: '',
      });
      setSelectedAtlas('');
      setSelectedModel('');
    }
  }, [atlases, models]);

  return (
    <Container>
      <Box width="500px" margin="auto" sx={{ marginTop: '1%', marginBottom: '1%' }}>
        <Stepper activeStep={activeStep}>
          {steps.map((labelText, index) => (
            <Step index={index}>
              <StepButton color="inherit" onClick={() => handleStep(index)}>
                {labelText}
              </StepButton>
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
              steps={steps}
              setSelectedAtlas={handleAtlasSelection}
              setSelectedModel={setSelectedModel}
              setActiveStep={handleStep}
              compatibleModels={selectedAtlas ? selectedAtlas.compatibleModels : []}
              atlases={atlases}
              models={models}
              demos={demoDatasets}
              selectedDataset={selectedDemoDataset}
              setSelectedDataset={setSelectedDemoDataset}
              datasetIsSelected={demoDatasetIsSelected}
              setDatasetIsSelected={setDemoDatasetIsSelected}
            />
          )
          : (
            <UploadFilePage
              path={path}
              selectedAtlas={selectedAtlas}
              selectedModel={selectedModel}
              setActiveStep={handleStep}
              demos={demoDatasets}
              selectedDataset={selectedDemoDataset}
              setSelectedDataset={setSelectedDemoDataset}
              datasetIsSelected={demoDatasetIsSelected}
              setDatasetIsSelected={setDemoDatasetIsSelected}
            />
          )
      }
    </Container>
  );
}

export default GeneMapperState;
