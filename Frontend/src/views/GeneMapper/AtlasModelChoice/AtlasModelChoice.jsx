import {
  Grid, Typography, Stack, Alert, Box, Tooltip,
} from '@mui/material';
import AtlasCardSelect from 'components/Cards/AtlasCardSelect';
import { DemoDatasetCard } from 'components/GeneMapper/TabCard';
import { ModelCardSelect } from 'components/Cards/ModelCardSelect';
import CustomButton from 'components/CustomButton';
import { colors } from 'shared/theme/colors';
import { useHistory } from 'react-router-dom';
import React, { useState, useEffect } from 'react';
import Clear from '@mui/icons-material/Clear';
import ArrowForwardIcon from '@mui/icons-material/ArrowForward';

function AtlasModelChoice({
  setActiveStep,
  selectedAtlas, setSelectedAtlas,
  selectedModel, setSelectedModel, path,
  compatibleModels, atlases, models,
}) {
  const [showWarning, setShowWarning] = useState(false);
  const history = useHistory();
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
      demoId: 123131,
      title: 'pbmc + totalvi',
      atlas: 'pbmc',
      model: 'totalvi',
      url: 'link to the dataset to fetch',
    },
  ];

  const [datasetIsSelected, setDatasetIsSelected] = useState(false);
  const [selectedDataset, setSelectedDataset] = useState(null);
  const [demoDatasets, setDemoDatasets] = useState(demoDatasetArr);

  // Make sure the choice here works

  const PBMCObj = {
    _id: '626ea3311d7d1a27de465b64',
    name: 'PBMC',
    previewPictureURL: 'https://storage.googleapis.com/jst-2021-bucket-static/pbmc.png',
    modalities: 'CITE-seq: CD3_TotalSeqB, CD4_TotalSeqB, CD8a_TotalSeqB, CD14_TotalSeqB, CD15_TotalSeqB, CD16_TotalSeqB, CD56_TotalSeqB, CD19_TotalSeqB, CD25_TotalSeqB, CD45RA_TotalSeqB, CD45RO_TotalSeqB, PD-1_TotalSeqB, TIGIT_TotalSeqB, CD127_TotalSeqB',
    numberOfCells: '11K',
    species: ['Human'],
    compatibleModels: ['totalVI'],
  };

  console.log(`The atlases are: ${JSON.stringify(atlases)}`);
  /* Handling the choice of the demo dataset */
  const handleDemoClick = (dataset) => {
    // only set selected dataset if not already selected or null
    if (!selectedDataset || selectedDataset.demoId !== dataset.demoId) {
      setSelectedDataset(dataset);
      setDatasetIsSelected(true);
      console.log(`The first element[0] ${JSON.stringify(atlases[0].name)}`);
      console.log(`The first element[1] ${JSON.stringify(atlases[1])}`);
      console.log(`The first element[2] ${JSON.stringify(atlases[2])}`);
      console.log(`The dataset atlas is: ${dataset.atlas}`);
      // find the objects corresponding to the atlas and model
      const atlasObj = atlases.filter(
        (a) => a.name.toLowerCase() === dataset.atlas.toLowerCase(),
      )[0];
      const modelObj = models.filter(
        (m) => m.name.toLowerCase() === dataset.model.toLowerCase(),
      )[0];

      setSelectedAtlas(atlasObj);
      setSelectedModel(modelObj);
    } else {
      setDatasetIsSelected(false);
      setSelectedDataset(null);
      // deselect the atlas and model
      setSelectedAtlas({});
      setSelectedModel({});
    }
  };

  /* fetch demo datasets */

  return (
    <div>
      {showWarning
        && (
          <Alert severity="error">
            Select an Atlas and a fitting Model before continuing
          </Alert>
        )}
      <Typography
        variant="h5"
        sx={{
          fontWeight: 'bold',
          pb: '1em',
          mt: '1.5em',
        }}
      >
        Pick an Atlas
      </Typography>

      <Grid container spacing={2} width="100%" overflow="auto" wrap="nowrap">
        {
          atlases && atlases.map((a) => (
            <Grid item height="320px">
              <AtlasCardSelect
                width="225px"
                height="97%"
                title={a.name}
                modalities={a.modalities}
                cellsInReference={a.numberOfCells}
                species={a.species}
                imgLink={a.previewPictureURL}
                selected={selectedAtlas.name === a.name}
                onSelect={setSelectedAtlas}
                atlasObject={a}
              />
            </Grid>
          ))
        }
      </Grid>

      <Typography variant="h5" sx={{ fontWeight: 'bold', pb: '1em' }} marginTop="32px">
        Pick a Model
      </Typography>

      <Grid container spacing={2} width="100%" overflow="auto" wrap="nowrap">
        {
          models && models.map((m) => (
            <Grid item height="225px">
              <ModelCardSelect
                width="225px"
                height="97%"
                title={m.name}
                description={m.description}
                selected={selectedModel.name === m.name}
                onSelect={setSelectedModel}
                modelObject={m}
                disabled={!compatibleModels
                  || !compatibleModels.map(
                    (m) => m.toLowerCase(),
                  ).includes(m.name.toLowerCase()) || compatibleModels.length === 0}
              />
            </Grid>
          ))
        }
      </Grid>

      {/* Demo datasets */}
      <Box sx={{ width: 400 }}>
        <Typography variant="h5" sx={{ fontSize: '18px', fontWeight: 'bold', pb: '1em' }} marginTop="32px">
          Or try out one of the available demo datasets
        </Typography>
      </Box>
      {/* Loop over demo datasets */}
      { demoDatasets.length !== 0 && demoDatasets.map((dataset) => (
        <DemoDatasetCard
          title={dataset.title}
          atlas={dataset.atlas}
          handleOnClick={() => handleDemoClick(dataset)}
          // selected={datasetIsSelected && dataset.demoId === selectedDataset.demoiId}
          selected={datasetIsSelected && dataset.demoId === selectedDataset.demoId}
        />
      ))}

      <Stack direction="row" justifyContent="space-between" sx={{ marginTop: '50px', marginBottom: '3em' }}>
        <CustomButton type="tertiary" onClick={() => history.push(`${path}`)}>
          <Clear />
          &nbsp; Cancel
        </CustomButton>
        <Tooltip title={!selectedAtlas || !selectedModel ? 'Select an Atlas and a fitting Model before continuing' : ''} placement="top">
          <Box
            onClick={!selectedAtlas || !selectedModel ? () => setShowWarning(true) : () => { }}
          >
            <CustomButton
              type="primary"
              disabled={!selectedAtlas || !selectedModel}
              onClick={() => setActiveStep(1)}
            >
              Confirm&nbsp;
              <ArrowForwardIcon />
            </CustomButton>
          </Box>
        </Tooltip>
      </Stack>
    </div>
  );
}

export default AtlasModelChoice;
