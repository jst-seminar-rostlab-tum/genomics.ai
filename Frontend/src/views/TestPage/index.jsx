import {Switch, Route} from 'react-router-dom'

import {Box} from '@mui/material'
import {TabGroup} from 'components/Tab'

import {colors} from 'shared/theme/colors'
import { Redirect } from 'react-router-dom'

const exampleTabsInfo=[
    {
        label: "tab 1",
        path: "/test/tab_1"
    },
    {
        label: "tab 2",
        path: "/test/tab_2"
    },
    {
        label: "tab 3",
        path: "/test/tab_3"
    }
]

export default function TestPage() {
    return (
        <>
            <Box sx={{width: "500px", height: "500px", bgcolor: "white", margin: "auto", border: 2}}>
                <TabGroup tabsInfo={exampleTabsInfo}/>
                <Switch>
                    {
                        exampleTabsInfo.map((tabInfo, index)=>(<Route path={tabInfo.path} render={()=><Box>tab {index+1}</Box>}/>))
                    }
                    <Redirect to={exampleTabsInfo[0].path} />
                </Switch>
            </Box>

            <Box sx={{width: "500px", height: "500px", bgcolor: colors.primary[800], margin: "auto", border: 2}}>
                <TabGroup tabsInfo={exampleTabsInfo} darkBackground/>
                <Switch>
                    {
                        exampleTabsInfo.map((tabInfo, index)=>(<Route path={tabInfo.path} render={()=><Box sx={{color: "white"}}>tab {index+1}</Box>}/>))
                    }
                    <Redirect to={exampleTabsInfo[0].path} />
                </Switch>
            </Box>
        </>
    )
}
