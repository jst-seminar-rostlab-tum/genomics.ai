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
      demoId: 1232112,
      title: 'demo dataset1',
      atlas: 'Pancreas',
      model: 'SCVI',
      url: 'link to the dataset to fetch',
    },
    {
      demoId: 123131,
      title: 'demo dataset2',
      atlas: 'PBMC',
      model: 'SCANVI',
      url: 'link to the dataset to fetch',
    },
    {
      demoId: 131,
      title: 'demo dataset1',
      atlas: 'Fetal Immune',
      model: 'TOTALVI',
      url: 'link to the dataset to fetch',
    },
  ];

  const [datasetIsSelected, setDatasetIsSelected] = useState(false);
  const [selectedDataset, setSelectedDataset] = useState(false);
  const [demoDatasets, setDemoDatasets] = useState(demoDatasetArr);
  const data = { name: 'hello, this is a data object', visibility: 'true' };
  { /* fetch the datasets */ }

  { /* Handling the choice of the demo dataset */ }
  const onClick = ({
    demoId, title, atlas, mode, url,
  }) => {
    alert(demoId);
    // set the dataset to the one with the matching id
  };

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
            <Grid item height="150px">
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
      { demoDatasetArr.length !== 0 && demoDatasetArr.map((dataset) => (
        <DemoDatasetCard
          title={dataset.title}
          atlas={dataset.atlas}
          handleOnClick={() => onClick(dataset)}
          selected={selectedDataset && dataset.demoId === selectedDataset.demoiId}
        />
      ))}

      {/* what to render as the demo datasets */}
      {/* ( demoDatasets.length !== 0 ?  */}

      {/* ( demoDatasets.length !== 0 ? demoDatasets.map((dataset) => <DemoDatasetCard
            title = {dataset}
            atlas = {dataset}
            handleOnClick = {() => onClick(dataset)}
            selected = {selectedDataset && dataset.demoId === selectedDataset.demoiId}
            />) : <Alert severity="info"> No existing datasets available. </Alert>;

        ) */}

      <Stack direction="row" justifyContent="space-between" sx={{ marginTop: '20px', marginBottom: '3em' }}>
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
