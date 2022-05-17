import React, { useState, useEffect, useCallback } from 'react';
import {
  useHistory, useLocation,
  useRouteMatch,
} from 'react-router-dom';
import { Box, Stack } from '@mui/material';
import { TabGroup } from 'components/Tab';
import Search from 'components/Search';
import Filter from 'components/ExplorePageComponents/Filter';
import NavBar from 'components/NavBar';
import Breadcrumb from 'components/Breadcrumb';
import LoginForm from 'components/LoginForm';
import RegistrationForm from 'components/RegistrationForm';

import ModelsService from 'shared/services/Models.service';
import AtlasService from 'shared/services/Atlas.service';
import AtlasesGrid from 'components/Grids/AtlasesGrid';
import ModelsGrid from 'components/Grids/ModelsGrid';
import Mapper from 'components/Mapper';
import { applyModelFilters, applyAtlasFilters } from 'shared/utils/filter';
import ExploreRoutes from 'components/ExplorePageComponents/ExploreRoutes';

const tmpObj = [
  {
    label: 'ATLASES',
    path: '/explore/atlases',
  },
  {
    label: 'MODELS',
    path: '/explore/models',
  },
];

const Explore = () => {
  const [value, setValue] = useState(0);
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

  useEffect(() => {
    AtlasService.getAtlases()
      .then((newAtlases) => setAtlases(newAtlases))
      .catch((err) => console.log(err));
    ModelsService.getModels()
      .then((newModels) => setModels(newModels))
      .catch((err) => console.log(err));
  }, []);

  const searchedKeywordChangeHandler = (value) => {
    updateQueryParams('keyword', value);
  };

  useEffect(() => {
    if (selectedAtlas || selectedModel) setMapperVisible(true);
    if (!selectedAtlas && !selectedModel) setMapperVisible(false);
  }, [selectedAtlas, selectedModel]);

  const tabMenu = () => (
    <Box height="50px">
      <TabGroup value={value} setValue={setValue} tabsInfo={tmpObj} />
      {value === 0 ? (
        <AtlasesGrid
          atlases={applyAtlasFilters(atlases, searchedKeyword, searchParams)}
          path={path}
          setSelectedAtlas={setSelectedAtlas}
          selectedAtlas={selectedAtlas}
        />
      ) : null }
      {value === 1 ? (
        <ModelsGrid
          models={applyModelFilters(models, searchedKeyword, searchParams)}
          path={path}
          setSelectedModel={setSelectedModel}
          selectedModel={selectedModel}
        />
      ) : null }
    </Box>

  );

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

  const tmp_elems = useLocation().pathname.split('/');
  const elems = tmp_elems.map((elem, index) => {
    if (index === 3) {
      if (tmp_elems[2] === 'atlases') return atlases.filter((x) => x._id === elem)[0] ? atlases.filter((x) => x._id === elem)[0].name : elem;
      if (tmp_elems[2] === 'models') return models.filter((x) => x._id === elem)[0] ? models.filter((x) => x._id === elem)[0].name : elem;
    }
    return elem;
  });

  const executeScroll = () => history.push({ pathname: '/', state: { contact_us: true } });

  return (
    <Box
      sx={{
        display: 'flex',
        flexDirection: 'column',
        '::-webkit-scrollbar': {
          display: 'none',
        },
        height: '100vh',
        overflow: 'hidden',
      }}
    >
      {isLoginFormVisible && (
        <LoginForm visible={isLoginFormVisible} onClose={onLoginFormClosed} />
      )}
      {isRegistrationFormVisible && (
        <RegistrationForm
          visible={isRegistrationFormVisible}
          onClose={onRegistrationFormClosed}
        />
      )}

      <Box>
        <NavBar
          position="relative"
          onLoginClicked={onLoginClicked}
          onSignUpClicked={onSignUpClicked}
          executeScroll={executeScroll}
        />
      </Box>

      <Stack
        direction="row"
        sx={{
          alignSelf: 'center', width: '60%', marginTop: '2%', justifyContent: 'space-between',
        }}
      >
        <Breadcrumb elems={elems} fontSize={1} actions={{ explore: () => setValue(0) }} />
        <Box sx={{ alignSelf: 'center', width: '40%', marginBlock: '2%' }}>
          <Search
            filterComponent={(
              <Filter
                searchParams={searchParams}
                updateQueryParams={updateQueryParams}
                path={path}
              />
            )}
            handleSearch={searchedKeywordChangeHandler}
            value={searchedKeyword}
            padding="0px"
          />
        </Box>
      </Stack>

      <Box
        sx={{
          display: 'flex',
          flexDirection: 'column',
          alignSelf: 'center',
          width: { xs: '90%', md: '60%' },
        }}
      >
        {/* /explore/atlases */}
        <ExploreRoutes atlases={atlases && tabMenu()} models={models && tabMenu()} path="/explore" />
      </Box>

      <Mapper
        mapperAtlas={selectedAtlas || null}
        mapperModel={selectedModel || null}
        setSelectedAtlas={setSelectedAtlas}
        setSelectedModel={setSelectedModel}
        open={mapperVisible}
        fabOnClick={() => setMapperVisible(!mapperVisible)}
      />
    </Box>
  );
};

export default Explore;
