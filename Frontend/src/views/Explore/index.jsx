import React, { useState, useEffect, useCallback } from 'react';
import {
  Switch, Route, Redirect, useHistory,
} from 'react-router-dom';
import { Box } from '@mui/material';
import Grid from '@mui/material/Grid';
import { TabGroup } from 'components/Tab';
import Search from 'components/Search';
import { Filter } from 'components/Filter/Filter';
import NavBar from 'components/NavBar';
import Breadcrumb from 'components/Breadcrumb';
import AtlasCard from 'components/Cards/AtlasCard';
import { ModelCard } from 'components/Cards/ModelCard';
import Mapper from 'components/Mapper';
import LoginForm from 'components/LoginForm';
import RegistrationForm from 'components/RegistrationForm';
import './Explore.css';

import DatasetCard from 'components/Cards/DatasetCard';
import ModelsService from 'shared/services/Models.service'
import AtlasService from 'shared/services/Atlas.service'

import { useLocation } from 'react-router-dom';
import { useAuth } from 'shared/context/authContext';

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
    label: "DATASETS",
    path: "/explore/datasets"
  }
]


const datasets = [
  {
    title: "dataset 1",
    category: "human"
  },
  {
    title: "dataset 2",
    category: "animal"
  }
]
const datasetsGrid = (
  <Box className="cardsContainer">
    <Grid container spacing={3}>
      {datasets.map((d, index) => (
        <Grid key={index} item xs={12} sm={6} md={4} lg={3}>
          <DatasetCard title={d.title} category={d.category} />
        </Grid>
      ))}
    </Grid>
  </Box>
)

const Explore = () => {

  const [value, setValue] = useState(0);
  const [selectedAtlas, setSelectedAtlas] = useState(null);
  const [selectedModel, setSelectedModel] = useState(null);
  const [mapperVisible, setMapperVisible] = useState(false);
  const [isLoginFormVisible, setLoginFormVisible] = useState(false);
  const [isRegistrationFormVisible, setRegistrationFormVisible] = useState(false);
  const [searchValue, setSearchValue] = useState('');
  const history = useHistory();

  const handleSearch = (e) => {
    setSearchValue(e);
    console.log(e);
  };

  const [atlases, setAtlases] = useState([])
  const [models, setModels] = useState([])

  const path = useLocation().pathname

  useEffect(() => {
    AtlasService.getAtlases()
      .then((atlases) => setAtlases(atlases))
      .catch((err) => console.log(err))

    ModelsService.getModels()
      .then((models) => setModels(models))
      .catch((err) => console.log(err))

  }, [])


  useEffect(() => {
    if(selectedAtlas || selectedModel) setMapperVisible(true)
    if(!selectedAtlas && !selectedModel) setMapperVisible(false)
  }, [selectedAtlas, selectedModel])

  const atlasesGrid = (atlases, path) => (
    <Box className='atlasContainer'>
      <Grid container spacing={3}>
        {
          atlases.map((atlas) => (
            <Grid key={atlas._id} item xs={12} sm={6} md={4} lg={3}>
              <AtlasCard onClick={() => setSelectedAtlas(atlas)} width="300px" height="500px" imgLink={atlas.previewPictureURL} species={atlas.species} modalities={atlas.modalities} title={atlas.name} learnMoreLink={`${path}/${atlas._id}`} />
            </Grid>
          ))
        }
      </Grid>
    </Box >
  )

  const modelsGrid = (models, path) => (
    <Box className='cardsContainer'>
      <Grid container spacing={3}>
        {
          models.map((model) => (
            <Grid key={model._id} item xs={12} sm={6} md={4} lg={3}>
              <ModelCard onClick={() => setSelectedModel(model)} title={`Model ${model.name}`} description={model.description} learnMoreLink={`${path}/${model._id}`} />
            </Grid>
          ))
        }
      </Grid>
    </Box>
  )
  const onLoginClicked = useCallback(() => {
    console.log("login")
    setRegistrationFormVisible(false)
    setLoginFormVisible(true)
  }, [setLoginFormVisible])

  const onSignUpClicked = useCallback(() => {
    console.log('register');
    setLoginFormVisible(false);
    setRegistrationFormVisible(true);
  }, [setRegistrationFormVisible]);

  const onLoginFormClosed = useCallback(() => {
    setLoginFormVisible(false);
  }, [setLoginFormVisible]);

  const onRegistrationFormClosed = useCallback(() => {
    setRegistrationFormVisible(false);
  }, [setRegistrationFormVisible]);

  const handleModelMapClick = (index) => {
    setSelectedModel(models.find(model => model._id === index));
    if (!mapperVisible) {
      setMapperVisible(true);
    }
    if (selectedAtlas === null) {
      history.push('/explore/atlases');
    }
  };

  return (
    <Box sx={{ display: 'flex', flexDirection: 'column', justifyContent: 'center' }}>

      {isLoginFormVisible
        && <LoginForm visible={isLoginFormVisible} onClose={onLoginFormClosed} />}
      {isRegistrationFormVisible
        && (
          <RegistrationForm
            visible={isRegistrationFormVisible}
            onClose={onRegistrationFormClosed}
          />
        )}

      <Box>
        <NavBar position="relative" onLoginClicked={onLoginClicked} onSignUpClicked={onSignUpClicked} />
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
          <Route path="/explore/atlases" render={() => atlasesGrid(atlases, path)} />
          <Route path="/explore/models" render={() => modelsGrid(models, path)} />
          <Route path="/explore/datasets" render={() => datasetsGrid} />
          <Redirect to="/explore/atlases" />
        </Switch>

      </Box>

      <Mapper
        mapperAtlas={selectedAtlas ? selectedAtlas.name : null}
        mapperModel={selectedModel ? selectedModel.name : null}
        setSelectedAtlas={setSelectedAtlas}
        setSelectedModel={setSelectedModel}
        open={mapperVisible}
        fabOnClick={() => setMapperVisible(!mapperVisible)}
      />
    </Box>
  );
}

export default Explore;
