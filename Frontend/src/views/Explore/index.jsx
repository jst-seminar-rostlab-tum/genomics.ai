import React from 'react'
import { Switch, Route, Redirect } from 'react-router-dom'
import { Box } from '@mui/material'
import { useLocation } from 'react-router-dom/cjs/react-router-dom.min'
import Grid from '@mui/material/Grid';

import { TabGroup } from 'components/Tab'
import Search from 'components/Search'
import Filter from 'components/Filter'
import NavBar from 'components/NavBar'
import Breadcrumb from "components/Breadcrumb"
import AtlasCard from 'components/Cards/AtlasCard'
import DatasetCard from 'components/Cards/DatasetCard'
import ModelCard from 'components/Cards/ModelCard'
import Mapper from "components/Mapper"

import './Explore.css'

const tmpObj = [
    {
        label: "ATLASES",
        path: "/explore/atlases"
    },
    {
        label: "MODELS",
        path: "/explore/models"
    },
    {
        label: "DATASETS",
        path: "/explore/datasets"
    }
]

const atlasesGrid = (
    <Box className='atlasContainer'>
        <Grid container spacing={2}>
            <Grid item xs={12} sm={6} md={4} lg={3}>
                <AtlasCard title='Human - PBMC' imgLink='https://3-h.de/wp-content/uploads/grey-gradient-background.jpeg' modalities=' RNA, ADT' cellsInReference='161,764' species='Human' />
            </Grid>
            <Grid item xs={12} sm={6} md={4} lg={3}>
                <AtlasCard title='Human - PBMC' imgLink='https://3-h.de/wp-content/uploads/grey-gradient-background.jpeg' modalities=' RNA, ADT' cellsInReference='161,764' species='Human' />
            </Grid>
            <Grid item xs={12} sm={6} md={4} lg={3}>
                <AtlasCard title='Human - PBMC' imgLink='https://3-h.de/wp-content/uploads/grey-gradient-background.jpeg' modalities=' RNA, ADT' cellsInReference='161,764' species='Human' />
            </Grid>
            <Grid item xs={12} sm={6} md={4} lg={3}>
                <AtlasCard title='Human - PBMC' imgLink='https://3-h.de/wp-content/uploads/grey-gradient-background.jpeg' modalities=' RNA, ADT' cellsInReference='161,764' species='Human' />
            </Grid>
            <Grid item xs={12} sm={6} md={4} lg={3}>
                <AtlasCard title='Human - PBMC' imgLink='https://3-h.de/wp-content/uploads/grey-gradient-background.jpeg' modalities=' RNA, ADT' cellsInReference='161,764' species='Human' />
            </Grid>
            <Grid item xs={12} sm={6} md={4} lg={3}>
                <AtlasCard title='Human - PBMC' imgLink='https://3-h.de/wp-content/uploads/grey-gradient-background.jpeg' modalities=' RNA, ADT' cellsInReference='161,764' species='Human' />
            </Grid>
            <Grid item xs={12} sm={6} md={4} lg={3}>
                <AtlasCard title='Human - PBMC' imgLink='https://3-h.de/wp-content/uploads/grey-gradient-background.jpeg' modalities=' RNA, ADT' cellsInReference='161,764' species='Human' />
            </Grid>
            <Grid item xs={12} sm={6} md={4} lg={3}>
                <AtlasCard title='Human - PBMC' imgLink='https://3-h.de/wp-content/uploads/grey-gradient-background.jpeg' modalities=' RNA, ADT' cellsInReference='161,764' species='Human' />
            </Grid>
        </Grid>
    </Box >
)

const modelsGrid = (
    <Box className='cardsContainer'>
        <Grid container spacing={2}>
            <Grid item xs={12} sm={6} md={4} lg={3}>

            </Grid>
            <Grid item xs={12} sm={6} md={4} lg={3}>

            </Grid>
            <Grid item xs={12} sm={6} md={4} lg={3}>

            </Grid>
        </Grid>
    </Box>
)

const datasetsGrid = (
    <Box className='cardsContainer'>
        <Grid container spacing={2}>
            <Grid item xs={12} sm={6} md={4} lg={3}>
                <DatasetCard />
            </Grid>
            <Grid item xs={12} sm={6} md={4} lg={3}>
                <DatasetCard />
            </Grid>
            <Grid item xs={12} sm={6} md={4} lg={3}>
                <DatasetCard />
            </Grid>
        </Grid>
    </Box>
)

const Explore = () => {

    const location = useLocation()

    return (
        <Box sx={{ display: 'flex', flexDirection: 'column', justifyContent: 'center' }} >
            <Box>
                <NavBar />
                <h1>NavBar goes here</h1>
            </Box>

            <Box sx={{ alignSelf: 'center', width: '65%', marginTop: '2%' }}>
                <Breadcrumb path={`${location.pathname}`} fontSize={1} />
            </Box>

            <Box sx={{ alignSelf: 'center', width: '65%', marginBlock: '2%' }}>
                <Search filterComponent={<Filter references={["test", "test"]} categories={["category1", "category2"]} />} handleSearch={(textRef) => console.log(textRef)} />
            </Box>

            <Box sx={{ display: 'flex', flexDirection: 'column', justifyContent: 'center', alignSelf: 'center', width: '80%' }}>
                {/* /explore/atlases */}
                <TabGroup tabsInfo={tmpObj} />

                <Switch>
                    <Route path="/explore/atlases" render={() => atlasesGrid} />
                    <Route path="/explore/models" render={() => modelsGrid} />
                    <Route path="/explore/datasets" render={() => datasetsGrid} />
                    <Redirect to="/explore/atlases" />
                </Switch>

            </Box>

            <Mapper />
        </Box>
    )
}

export default Explore