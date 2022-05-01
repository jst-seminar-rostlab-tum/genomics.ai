import React, { useState } from 'react';
import {
  Switch, Route, Redirect, useHistory,
} from 'react-router-dom';
import { Box } from '@mui/material';
import Grid from '@mui/material/Grid';

import { TabGroup } from 'components/Tab';
import Search from 'components/Search';

import Filter from 'components/Filter/old';
import NavBar from 'components/NavBar';
import Breadcrumb from 'components/Breadcrumb';
import AtlasCard from 'components/Cards/AtlasCard';
import DatasetCard from 'components/Cards/DatasetCard';
import { ModelCard } from 'components/Cards/ModelCard';
import Mapper from 'components/Mapper';

import './Explore.css';
import { FunctionsRounded } from '@mui/icons-material';

const tmpObj = [
  {
    label: 'ATLASES',
    path: '/explore/atlases',
  },
  {
    label: 'MODELS',
    path: '/explore/models',
  },
  {
    label: 'DATASETS',
    path: '/explore/datasets',
  },
];

const atlasTitles = ['Human - PBMC0', 'Human - PBMC1', 'Human - PBMC2', 'Human - PBMC3', 'Human - PBMC4', 'Human - PBMC5', 'Human - PBMC6', 'Human - PBMC7'];
const atlasImgLinks = ['https://3-h.de/wp-content/uploads/grey-gradient-background.jpeg',
  'https://3-h.de/wp-content/uploads/grey-gradient-background.jpeg',
  'https://3-h.de/wp-content/uploads/grey-gradient-background.jpeg',
  'https://3-h.de/wp-content/uploads/grey-gradient-background.jpeg',
  'https://3-h.de/wp-content/uploads/grey-gradient-background.jpeg',
  'https://3-h.de/wp-content/uploads/grey-gradient-background.jpeg',
  'https://3-h.de/wp-content/uploads/grey-gradient-background.jpeg',
  'https://3-h.de/wp-content/uploads/grey-gradient-background.jpeg'];
const atlasModalities = ['RNA, ADT', 'RNA, ADT', 'RNA, ADT', 'RNA, ADT', 'RNA, ADT', 'RNA, ADT', 'RNA, ADT', 'RNA, ADT'];
const atlasCellsInReference = ['161,764', '161,764', '161,764', '161,764', '161,764', '161,764', '161,764', '161,764'];
const atlasSpecies = ['Human', 'Human', 'Human', 'Human', 'Human', 'Human', 'Human', 'Human'];

const modelTitles = ['Model 1', 'Model 2', 'Model 3', 'Model 4'];
const modelDescriptions = ['This is a short description', 'This is a short description', 'This is a short description', 'This is a short description', 'This is a short description'];

const datasetTitles = ['title1', 'title2', 'title3'];
const datasetCategories = ['category1', 'category1', 'category1'];
const datasetsGrid = (
  <Box className="cardsContainer">
    <Grid container spacing={3}>
      {datasetTitles.map((title, index) => (
        <Grid item xs={12} sm={6} md={4} lg={3}>
          <DatasetCard title={title} category={datasetCategories[index]} />
        </Grid>
      ))}
    </Grid>
  </Box>
);

const Explore = () => {
  const [value, setValue] = useState(0);
  const [selectedAtlas, setSelectedAtlas] = useState(-1);
  const [selectedModel, setSelectedModel] = useState(-1);
  const [mapperVisible, setMapperVisible] = useState(false);
  const [searchValue, setSearchValue] = useState('');
  const history = useHistory();

  const handleSearch = (e) => {
    setSearchValue(e);
    console.log(e);
  };

  const handleAtlasMapClick = (index) => {
    setSelectedAtlas(index);
    if (!mapperVisible) {
      setMapperVisible(true);
    }
    if (selectedModel === -1) {
      history.push('/explore/models');
    }
  };

  const handleModelMapClick = (index) => {
    setSelectedModel(index);
    if (!mapperVisible) {
      setMapperVisible(true);
    }
    if (selectedAtlas === -1) {
      history.push('/explore/atlases');
    }
  };

  const atlasesGrid = (
    <Box className="atlasContainer">
      <Grid container spacing={3}>
        {atlasTitles.map((title, index) => (
          <Grid item xs={12} sm={6} md={4} lg={3}>
            <AtlasCard width="300px" height="500px" index={index} title={title} imgLink={atlasImgLinks[index]} modalities={atlasModalities[index]} cellsInReference={atlasCellsInReference[index]} species={atlasSpecies[index]} onClick={handleAtlasMapClick} />
          </Grid>
        ))}
      </Grid>
    </Box>
  );

  const modelsGrid = (
    <Box className="cardsContainer">
      <Grid container spacing={3}>
        {modelTitles.map((title, index) => (
          <Grid item xs={12} sm={6} md={4} lg={3}>
            <ModelCard
              title={title}
              description={modelDescriptions[index]}
              onClick={handleModelMapClick}
              index={index}
            />
          </Grid>
        ))}
      </Grid>
    </Box>
  );

  return (
    <Box sx={{ display: 'flex', flexDirection: 'column', justifyContent: 'center' }}>

      <Box>
        <NavBar />
        <h1>NavBar goes here</h1>
      </Box>

      <Box sx={{ alignSelf: 'center', width: '65%', marginTop: '2%' }}>
        <Breadcrumb fontSize={1} actions={{ explore: () => setValue(0) }} />
      </Box>

      <Box sx={{ alignSelf: 'center', width: '65%', marginBlock: '2%' }}>
        <Search filterComponent={<Filter references={['test', 'test']} categories={['category1', 'category2']} />} handleSearch={handleSearch} value={searchValue} />
      </Box>

      <Box sx={{
        display: 'flex', flexDirection: 'column', justifyContent: 'center', alignSelf: 'center', width: '80%',
      }}
      >
        {/* /explore/atlases */}
        <TabGroup value={value} setValue={setValue} tabsInfo={tmpObj} />
        <Switch>
          <Route path="/explore/atlases" render={() => atlasesGrid} />
          <Route path="/explore/models" render={() => modelsGrid} />
          <Route path="/explore/datasets" render={() => datasetsGrid} />
          <Redirect to="/explore/atlases" />
        </Switch>

      </Box>

      <Mapper
        mapperAtlas={atlasTitles[selectedAtlas] ?? null}
        mapperModel={modelTitles[selectedModel] ?? null}
        setSelectedAtlas={setSelectedAtlas}
        setSelectedModel={setSelectedModel}
        open={mapperVisible}
        fabOnClick={() => setMapperVisible(!mapperVisible)}
      />
    </Box>
  );
};

export default Explore;
