import React, {
  useState, useEffect, useCallback, useContext,
} from 'react';
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
import Footer from 'components/Footer';

import ModelsService from 'shared/services/Models.service';
import AtlasService from 'shared/services/Atlas.service';
import AtlasesGrid from 'components/Grids/AtlasesGrid';
import ModelsGrid from 'components/Grids/ModelsGrid';
import Mapper from 'components/Mapper';
import { applyModelFilters, applyAtlasFilters } from 'shared/utils/filter';
import ExploreRoutes from 'components/ExplorePageComponents/ExploreRoutes';
import { useAuth } from 'shared/context/authContext';
import { LoginContext } from 'shared/context/loginContext';
import PasswordForgetForm from 'components/PasswordForgetForm';

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
  const { search, pathname } = useLocation();
  const searchParams = new URLSearchParams(search);
  const searchedKeyword = searchParams.get('keyword') || '';
  const { path } = useRouteMatch();
  const [atlases, setAtlases] = useState([]);
  const [models, setModels] = useState([]);
  const [user, setUser] = useAuth();
  const history = useHistory();

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

  const handleAtlasSelection = (newAtlas) => {
    setSelectedAtlas(newAtlas);
    if (!selectedModel) {
      history.push(`${path}/models`);
      setValue(1);
    }
  };

  const handleModelSelection = (newModel) => {
    setSelectedModel(newModel);
    if (!selectedAtlas) {
      history.push(`${path}/atlases`);
      setValue(0);
    }
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

  const onValueChange = (newValue) => {
    setValue(newValue);
    searchedKeywordChangeHandler('');
  };

  const handleMap = () => {
    if (user) {
      history.push({ pathname: `/sequencer/genemapper/create/${selectedAtlas._id}/${selectedModel._id}` });
      setMapperVisible(false);
    }
  };

  useEffect(() => {
    if (selectedAtlas || selectedModel) setMapperVisible(true);
    if (!selectedAtlas && !selectedModel) setMapperVisible(false);
  }, [selectedAtlas, selectedModel]);

  const tabMenu = () => (
    <Box height="50px">

      <TabGroup value={value} onValueChange={onValueChange} tabsInfo={tmpObj} />
      {value === 0 ? (
        <AtlasesGrid
          atlases={applyAtlasFilters(atlases, searchedKeyword, searchParams, selectedModel)}
          path={path}
          handleAtlasSelection={handleAtlasSelection}
          selectedAtlas={selectedAtlas}
          selectedModel={selectedModel}
        />
      ) : null}
      {value === 1 ? (
        <ModelsGrid
          models={applyModelFilters(models, searchedKeyword, searchParams, selectedAtlas)}
          path={path}
          handleModelSelection={handleModelSelection}
          selectedModel={selectedModel}
          compatibleModels={selectedAtlas && selectedAtlas.compatibleModels}
        />
      ) : null}
    </Box>

  );

  const context = useContext(LoginContext);

  const onLoginClicked = () => {
    context.switchRegister(false);
    context.switchLogin(true);
  };

  const onSignUpClicked = () => {
    context.switchLogin(false);
    context.switchRegister(true);
  };

  const tmp_elems = pathname.split('/');
  const elems = tmp_elems.map((elem, index) => {
    if (index === 3) {
      if (tmp_elems[2] === 'atlases') return atlases.filter((x) => x._id === elem)[0] ? atlases.filter((x) => x._id === elem)[0].name : elem;
      if (tmp_elems[2] === 'models') return models.filter((x) => x._id === elem)[0] ? models.filter((x) => x._id === elem)[0].name : elem;
    }
    return elem;
  });

  const executeScroll = () => (user ? history.push({ pathname: '/sequencer/help' }) : history.push({ pathname: '/', state: { contact_us: true } }));

  return (
    <>
      <Box
        sx={{
          display: 'flex',
          flexDirection: 'column',
          '::-webkit-scrollbar': {
            display: 'none',
          },
          height: '100vh',
          overflow: 'auto',
        }}
      >
        {context.loginVisible && <LoginForm />}
        {context.registerVisible && <RegistrationForm />}
        {context.forgetVisible && <PasswordForgetForm />}

        <Box>
          <NavBar
            position="relative"
            onLoginClicked={onLoginClicked}
            onSignUpClicked={onSignUpClicked}
            executeScroll={executeScroll}
          />
        </Box>

        <Stack
          direction={{ xs: 'column', sm: 'column', md: 'row' }}
          sx={{
            alignSelf: 'center', width: { sx: '90%', md: '60%' }, marginTop: '2%', justifyContent: 'space-between',
          }}
        >
          <Breadcrumb elems={elems} fontSize={1} actions={{ explore: () => setValue(0) }} />
          <Box sx={{ alignSelf: 'center', width: { xs: '100%', md: '60%' }, marginBlock: '2%' }}>
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
              visible={pathname.split('/').slice(-1).includes('atlases') || pathname.split('/').slice(-1).includes('models')}
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
          <ExploreRoutes atlases={atlases && tabMenu()} models={models && tabMenu()} path="/explore" handleSelectAtlases={handleAtlasSelection} handleSelectModels={handleModelSelection} />
        </Box>

        <Mapper
          mapperAtlas={selectedAtlas ? selectedAtlas.name : null}
          mapperModel={selectedModel ? selectedModel.name : null}
          handleAtlasSelection={handleAtlasSelection}
          handleModelSelection={handleModelSelection}
          open={mapperVisible}
          fabOnClick={() => setMapperVisible(!mapperVisible)}
          handleMap={handleMap}
        />
      </Box>
      <Footer />
    </>
  );
};

export default Explore;
