import React, { useState, useEffect, useCallback } from 'react';
import {
  Switch, Route, Redirect, useHistory, useLocation, useParams,
  useRouteMatch,
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
import ModelsService from 'shared/services/Models.service';
import AtlasService from 'shared/services/Atlas.service';

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

const datasets = [
  {
    title: 'dataset 1',
    category: 'human',
  },
  {
    title: 'dataset 2',
    category: 'animal',
  },
];
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
);

const Explore = () => {
  const [tabValue, setTabValue] = useState(0);
  const [selectedAtlas, setSelectedAtlas] = useState(null);
  const [selectedModel, setSelectedModel] = useState(null);
  const [mapperVisible, setMapperVisible] = useState(false);
  const [isLoginFormVisible, setLoginFormVisible] = useState(false);
  const [isRegistrationFormVisible, setRegistrationFormVisible] = useState(false);
  const { search } = useLocation();
  const searchParams = new URLSearchParams(search);
  const searchedKeyword = searchParams.get('keyword') || '';
  const { path } = useRouteMatch();
  const history = useHistory();

  const [atlases, setAtlases] = useState([]);
  const [models, setModels] = useState([]);

  // function to update the state in the URL
  const updateQueryParams = (param, value) => {
    const params = new URLSearchParams(history.location.search);
    if (value) {
      params.set(param, value);
    } else {
      params.delete(param);
    }

    history.push({
      pathname: history.location.pathname,
      search: params.toString(),
    });
  };

  const searchedKeywordChangeHandler = (value) => {
    updateQueryParams('keyword', value);
  };

  useEffect(() => {
    AtlasService.getAtlases()
      .then((newAtlases) => setAtlases(newAtlases))
      .catch((err) => console.log(err));

    ModelsService.getModels()
      .then((newModels) => setModels(newModels))
      .catch((err) => console.log(err));
  }, []);

  useEffect(() => {
    if (selectedAtlas || selectedModel) setMapperVisible(true);
    if (!selectedAtlas && !selectedModel) setMapperVisible(false);
  }, [selectedAtlas, selectedModel]);

  function atlasesGrid() {
    const searchedAtlases = atlases.filter(
      (item) => item.name.toLowerCase().includes(searchedKeyword),
    );
    return (
      <Box className="atlasContainer" sx={{ height: '70vh' }}>
        <Grid container spacing={3}>
          {
            searchedAtlases.map((atlas) => (
              <Grid key={atlas._id} item xs={12} sm={6} md={4} lg={3}>
                <AtlasCard onClick={() => setSelectedAtlas(atlas)} imgLink={atlas.previewPictureURL} species={atlas.species} modalities={atlas.modalities} title={atlas.name} learnMoreLink={`${path}/${atlas._id}`} />
              </Grid>
            ))
          }
        </Grid>
      </Box>
    );
  }

  function modelsGrid() {
    const searchedModels = models.filter(
      (item) => item.name.toLowerCase().includes(searchedKeyword),
    );
    return (
      <Box className="cardsContainer" sx={{ height: '100%' }}>
        <Grid container spacing={3}>
          {
            searchedModels.map((model) => (
              <Grid key={model._id} item xs={12} sm={6} md={4} lg={3}>
                <ModelCard onClick={() => setSelectedModel(model)} title={`Model ${model.name}`} description={model.description} learnMoreLink={`${path}/${model._id}`} />
              </Grid>
            ))
          }
        </Grid>
      </Box>
    );
  }

  const onLoginClicked = useCallback(() => {
    console.log('login');
    setRegistrationFormVisible(false);
    setLoginFormVisible(true);
  }, [setLoginFormVisible]);

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
        <Breadcrumb fontSize={1} actions={{ explore: () => setTabValue(0) }} />
      </Box>

      <Box sx={{ alignSelf: 'center', width: '65%', marginBlock: '2%' }}>
        <Search filterComponent={<Filter references={['test', 'test']} categories={['category1', 'category2']} />} handleSearch={searchedKeywordChangeHandler} value={searchedKeyword} />
      </Box>

      <Box sx={{
        display: 'flex', flexDirection: 'column', alignSelf: 'center', width: '80%',
      }}
      >
        {/* /explore/atlases */}
        <TabGroup value={tabValue} setValue={setTabValue} tabsInfo={tmpObj} />
        <Switch>
          <Route path="/explore/atlases" render={() => atlasesGrid()} />
          <Route path="/explore/models" render={() => modelsGrid()} />
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
};

export default Explore;
